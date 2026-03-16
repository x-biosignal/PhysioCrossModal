# Signal alignment functions for PhysioCrossModal
#
# Functions for resampling and aligning multiple PhysioExperiment objects
# to a common sampling rate, and merging channels from different modalities.

# ---- Resampling internals ---------------------------------------------------

#' Resample a matrix to a new number of rows
#'
#' @param data Numeric matrix (time x channels).
#' @param n_old Integer number of original time points.
#' @param n_new Integer number of target time points.
#' @param method Character resampling method.
#' @return Resampled numeric matrix.
#' @noRd
.resample_matrix <- function(data, n_old, n_new, method) {
  duration <- 1.0  # normalised; only ratios matter
  t_old <- seq(0, duration, length.out = n_old)
  t_new <- seq(0, duration, length.out = n_new)

  n_channels <- ncol(data)
  result <- matrix(NA_real_, nrow = n_new, ncol = n_channels)

  for (ch in seq_len(n_channels)) {
    vec <- data[, ch]
    if (method == "linear") {
      result[, ch] <- stats::approx(t_old, vec, xout = t_new,
                                     method = "linear", rule = 2)$y
    } else if (method == "spline") {
      result[, ch] <- stats::spline(t_old, vec, xout = t_new,
                                     method = "natural")$y
    } else if (method == "fft") {
      result[, ch] <- .resample_fft(vec, n_new)
    }
  }
  result
}

#' FFT-based resampling for a single vector
#'
#' @param x Numeric vector.
#' @param n_new Target length.
#' @return Resampled numeric vector of length \code{n_new}.
#' @noRd
.resample_fft <- function(x, n_new) {
  n_old <- length(x)
  if (n_new == n_old) return(x)

  fft_x <- stats::fft(x)

  if (n_new > n_old) {
    # Upsample: zero-pad in frequency domain
    half <- floor(n_old / 2)
    fft_new <- complex(n_new)
    # Copy positive frequencies (DC through just below Nyquist)
    fft_new[seq_len(half)] <- fft_x[seq_len(half)]
    # Copy negative frequencies
    fft_new[(n_new - half + 1):n_new] <- fft_x[(n_old - half + 1):n_old]
    # Handle Nyquist bin for even-length signals: split between pos and neg
    if (n_old %% 2 == 0) {
      nyquist_val <- fft_x[half + 1] / 2
      fft_new[half + 1] <- nyquist_val
      fft_new[n_new - half + 1] <- nyquist_val
    } else {
      fft_new[half + 1] <- fft_x[half + 1]
    }
  } else {
    # Downsample: truncate in frequency domain
    half <- floor(n_new / 2)
    fft_new <- complex(n_new)
    # Copy positive frequencies (DC through half)
    fft_new[seq_len(half + 1)] <- fft_x[seq_len(half + 1)]
    # Copy negative frequencies from end of original
    n_neg <- n_new - half - 1L
    if (n_neg > 0L) {
      fft_new[(n_new - n_neg + 1):n_new] <- fft_x[(n_old - n_neg + 1):n_old]
    }
  }

  # R's fft(inverse=TRUE) does NOT divide by N, so divide by n_new
  # then scale by n_new/n_old for amplitude preservation => / n_old
  Re(stats::fft(fft_new, inverse = TRUE)) / n_old
}

# ---- alignToRate ------------------------------------------------------------

#' Resample a PhysioExperiment to a target sampling rate
#'
#' Resamples all assays in a \code{PhysioExperiment} object so that the
#' resulting data correspond to the specified target sampling rate. If the
#' \pkg{PhysioPreprocess} package is available and \code{method = "linear"},
#' its \code{resample()} function is used; otherwise an internal
#' interpolation is applied.
#'
#' @param x A \code{PhysioExperiment} object.
#' @param target_rate Numeric scalar: desired sampling rate in Hz.
#' @param method Character: resampling method. One of \code{"linear"}
#'   (default), \code{"spline"}, or \code{"fft"}.
#' @return A new \code{PhysioExperiment} with resampled data and updated
#'   \code{samplingRate}. The assay name is preserved from the original.
#' @seealso [alignSignals()], [mergePhysio()],
#'   [MultiPhysioExperiment()]
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(1000 * 4), nrow = 1000, ncol = 4)
#' pe <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   colData = S4Vectors::DataFrame(
#'     label = paste0("Ch", 1:4),
#'     type = rep("EEG", 4)
#'   ),
#'   samplingRate = 1000
#' )
#' pe500 <- alignToRate(pe, target_rate = 500)
#' samplingRate(pe500)  # 500
alignToRate <- function(x, target_rate,
                        method = c("linear", "spline", "fft")) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate required for resampling", call. = FALSE)
  }
  if (!is.numeric(target_rate) || length(target_rate) != 1 ||
      target_rate <= 0) {
    stop("target_rate must be a single positive number", call. = FALSE)
  }


  # No resampling needed if rates already match
  if (abs(sr - target_rate) < 1e-6) {
    return(x)
  }

  # Try PhysioPreprocess::resample() when available and method is compatible
  if (method == "linear" &&
      requireNamespace("PhysioPreprocess", quietly = TRUE)) {
    return(PhysioPreprocess::resample(x, target_rate, method = "linear",
                                       output_assay = defaultAssay(x)))
  }

  # Internal resampling -------------------------------------------------------
  assay_names <- SummarizedExperiment::assayNames(x)
  if (length(assay_names) > 1) {
    warning("Only the first assay ('", assay_names[1],
            "') is resampled. Additional assays are dropped.")
  }
  data <- SummarizedExperiment::assay(x, 1L)
  dims <- dim(data)
  n_old <- dims[1]

  duration <- n_old / sr
  n_new <- as.integer(round(duration * target_rate))

  # Anti-aliasing: apply lowpass filter before downsampling
  if (target_rate < sr) {
    aa_cutoff <- 0.9 * target_rate / 2
    if (length(dims) == 2) {
      for (ch in seq_len(dims[2])) {
        data[, ch] <- .lowpass_filter(data[, ch], cutoff = aa_cutoff, sr = sr)
      }
    } else if (length(dims) == 3) {
      for (s in seq_len(dims[3])) {
        for (ch in seq_len(dims[2])) {
          data[, ch, s] <- .lowpass_filter(data[, ch, s], cutoff = aa_cutoff,
                                            sr = sr)
        }
      }
    }
  }

  if (length(dims) == 2) {
    resampled <- .resample_matrix(data, n_old, n_new, method)
  } else if (length(dims) == 3) {
    resampled <- array(NA_real_, dim = c(n_new, dims[2], dims[3]))
    for (s in seq_len(dims[3])) {
      resampled[, , s] <- .resample_matrix(data[, , s], n_old, n_new, method)
    }
  } else {
    stop("Only 2D or 3D assays are supported", call. = FALSE)
  }

  new_row_data <- S4Vectors::DataFrame(time_idx = seq_len(n_new))

  new_pe <- PhysioExperiment(
    assays = S4Vectors::SimpleList(resampled),
    rowData = new_row_data,
    colData = SummarizedExperiment::colData(x),
    metadata = S4Vectors::metadata(x),
    samplingRate = target_rate
  )

  names(SummarizedExperiment::assays(new_pe)) <- assay_names[1]
  new_pe
}

# ---- alignSignals -----------------------------------------------------------

#' Align multiple PhysioExperiment objects to a common sampling rate
#'
#' Takes two or more named \code{PhysioExperiment} objects and resamples them
#' so they all share a single sampling rate, then wraps the result in a
#' \code{\link{MultiPhysioExperiment}}.
#'
#' Three alignment strategies are available:
#' \describe{
#'   \item{\code{"lowest_rate"}}{(default) Resample all signals to the lowest
#'     sampling rate present.}
#'   \item{\code{"common_rate"}}{Resample all signals to the highest sampling
#'     rate present.}
#'   \item{\code{"resample"}}{Resample all signals to the rate given by
#'     \code{target_rate} (which must be supplied).}
#' }
#'
#' If every input already shares the same sampling rate, no resampling is
#' performed regardless of the chosen method.
#'
#' @param ... Named \code{PhysioExperiment} objects.
#' @param method Character: one of \code{"lowest_rate"} (default),
#'   \code{"common_rate"}, or \code{"resample"}.
#' @param target_rate Numeric scalar: required when \code{method = "resample"}.
#' @return A \code{\link{MultiPhysioExperiment}} containing the (possibly
#'   resampled) input objects.
#' @seealso [alignToRate()], [mergePhysio()],
#'   [MultiPhysioExperiment()], [couplingAnalysis()]
#' @export
#' @examples
#' eeg <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(500 * 2), nrow = 500)),
#'   samplingRate = 500
#' )
#' emg <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(1000 * 2), nrow = 1000)),
#'   samplingRate = 1000
#' )
#' mpe <- alignSignals(EEG = eeg, EMG = emg, method = "lowest_rate")
alignSignals <- function(...,
                         method = c("lowest_rate", "common_rate", "resample"),
                         target_rate = NULL) {
  method <- match.arg(method)
  pes <- list(...)

  if (length(pes) == 0) {
    stop("At least one PhysioExperiment must be provided", call. = FALSE)
  }
  if (is.null(names(pes)) || any(names(pes) == "")) {
    stop("All PhysioExperiment arguments must be named", call. = FALSE)
  }
  if (!all(vapply(pes, inherits, logical(1), "PhysioExperiment"))) {
    stop("All arguments in ... must be PhysioExperiment objects", call. = FALSE)
  }

  rates <- vapply(pes, samplingRate, numeric(1))

  # Determine target rate
  if (method == "resample") {
    if (is.null(target_rate)) {
      stop("target_rate is required when method = 'resample'", call. = FALSE)
    }
  } else if (method == "common_rate") {
    target_rate <- max(rates)
  } else {
    # lowest_rate (default)
    target_rate <- min(rates)
  }

  # If all rates already match the target, skip resampling
  all_match <- all(abs(rates - target_rate) < 1e-6)
  if (!all_match) {
    pes <- lapply(pes, function(pe) {
      if (abs(samplingRate(pe) - target_rate) < 1e-6) {
        pe
      } else {
        alignToRate(pe, target_rate)
      }
    })
  }

  do.call(MultiPhysioExperiment, pes)
}

# ---- mergePhysio ------------------------------------------------------------

#' Merge two PhysioExperiment objects by combining channels
#'
#' Horizontally concatenates the channels (columns) of two
#' \code{PhysioExperiment} objects that share the **same** sampling rate.
#' Channel labels are prefixed to avoid name collisions.
#'
#' Both objects must have the same number of time points (rows) and the same
#' sampling rate. Only the first assay of each object is merged.
#'
#' @param x A \code{PhysioExperiment} object.
#' @param y A \code{PhysioExperiment} object.
#' @param prefix Character vector of length 2: prefixes added to channel
#'   labels from \code{x} and \code{y} respectively (default
#'   \code{c("x_", "y_")}).
#' @return A single \code{PhysioExperiment} with the columns of both inputs.
#' @seealso [alignToRate()], [alignSignals()],
#'   [MultiPhysioExperiment()]
#' @export
#' @examples
#' pe1 <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(100 * 4), nrow = 100)),
#'   colData = S4Vectors::DataFrame(label = paste0("A", 1:4)),
#'   samplingRate = 500
#' )
#' pe2 <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(100 * 4), nrow = 100)),
#'   colData = S4Vectors::DataFrame(label = paste0("B", 1:4)),
#'   samplingRate = 500
#' )
#' merged <- mergePhysio(pe1, pe2)
mergePhysio <- function(x, y, prefix = c("x_", "y_")) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(inherits(y, "PhysioExperiment"))

  if (length(prefix) != 2 || !is.character(prefix)) {
    stop("prefix must be a character vector of length 2", call. = FALSE)
  }

  sr_x <- samplingRate(x)
  sr_y <- samplingRate(y)
  if (abs(sr_x - sr_y) > 1e-6) {
    stop(
      sprintf(
        "Sampling rates must match for mergePhysio (got %.1f and %.1f). ",
        sr_x, sr_y
      ),
      "Use alignToRate() first to resample.",
      call. = FALSE
    )
  }

  data_x <- SummarizedExperiment::assay(x, 1L)
  data_y <- SummarizedExperiment::assay(y, 1L)

  if (nrow(data_x) != nrow(data_y)) {
    stop(
      sprintf(
        "Number of time points must match (got %d and %d)",
        nrow(data_x), nrow(data_y)
      ),
      call. = FALSE
    )
  }

  merged_data <- cbind(data_x, data_y)

  # Build merged colData with prefixed labels

  cd_x <- SummarizedExperiment::colData(x)
  cd_y <- SummarizedExperiment::colData(y)

  # Prefix the label column if it exists; otherwise create one
  if ("label" %in% colnames(cd_x)) {
    cd_x$label <- paste0(prefix[1], cd_x$label)
  } else {
    cd_x$label <- paste0(prefix[1], "ch", seq_len(ncol(data_x)))
  }
  if ("label" %in% colnames(cd_y)) {
    cd_y$label <- paste0(prefix[2], cd_y$label)
  } else {
    cd_y$label <- paste0(prefix[2], "ch", seq_len(ncol(data_y)))
  }

  # Harmonise colData columns before rbind
  all_cols <- union(colnames(cd_x), colnames(cd_y))
  for (col in all_cols) {
    if (!col %in% colnames(cd_x)) {
      cd_x[[col]] <- rep(NA, nrow(cd_x))
    }
    if (!col %in% colnames(cd_y)) {
      cd_y[[col]] <- rep(NA, nrow(cd_y))
    }
  }
  merged_cd <- rbind(cd_x[, all_cols, drop = FALSE], cd_y[, all_cols, drop = FALSE])

  # Combine metadata
  meta_x <- S4Vectors::metadata(x)
  meta_y <- S4Vectors::metadata(y)
  merged_meta <- c(meta_x, meta_y)

  assay_name_x <- SummarizedExperiment::assayNames(x)[1]
  assay_name <- if (!is.na(assay_name_x)) assay_name_x else "raw"

  new_pe <- PhysioExperiment(
    assays = S4Vectors::SimpleList(merged_data),
    colData = merged_cd,
    metadata = merged_meta,
    samplingRate = sr_x
  )

  names(SummarizedExperiment::assays(new_pe)) <- assay_name
  new_pe
}
