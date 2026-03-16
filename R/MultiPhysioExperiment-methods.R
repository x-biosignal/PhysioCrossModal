#' Additional methods for MultiPhysioExperiment
#'
#' Provides subsetting (`[`), `length()`, and `names()` methods for
#' \code{\link{MultiPhysioExperiment}} objects.

# ---- length method -------------------------------------------------------

#' Length of a MultiPhysioExperiment
#'
#' Returns the number of modalities (same as \code{\link{nModalities}}).
#'
#' @param x A \code{MultiPhysioExperiment} object.
#' @return Integer: number of modalities.
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   samplingRate = 250
#' )
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg)
#' length(mpe)  # 1
setMethod("length", "MultiPhysioExperiment", function(x) {
  nModalities(x)
})

# ---- names method --------------------------------------------------------

#' Names of a MultiPhysioExperiment
#'
#' Returns modality names (same as \code{\link{modalities}}).
#'
#' @param x A \code{MultiPhysioExperiment} object.
#' @return Character vector of modality names.
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
#' emg_data <- matrix(rnorm(1000 * 2), nrow = 1000, ncol = 2)
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   samplingRate = 250
#' )
#' pe_emg <- PhysioExperiment(
#'   assays = list(raw = emg_data),
#'   samplingRate = 1000
#' )
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)
#' names(mpe)  # c("EEG", "EMG")
setMethod("names", "MultiPhysioExperiment", function(x) {
  modalities(x)
})

# ---- [ subsetting method -------------------------------------------------

#' Subset a MultiPhysioExperiment by time range and/or modality
#'
#' Provides flexible subsetting of \code{MultiPhysioExperiment} objects:
#' \itemize{
#'   \item \code{mpe[, "eeg"]} -- select modalities by name, returns a new MPE
#'   \item \code{mpe[c(1.0, 5.0), ]} -- select a time window in seconds across
#'     all modalities, accounting for each modality's sampling rate
#'   \item \code{mpe[c(1.0, 5.0), "eeg"]} -- both time and modality subsetting
#' }
#'
#' When \code{i} is a numeric vector of length 2, it is interpreted as a time
#' range \code{[tmin, tmax]} in seconds.
#' For each modality, the appropriate sample indices are computed from its
#' sampling rate:
#' \code{start_idx = floor(tmin * sr) + 1}, \code{end_idx = floor(tmax * sr) + 1},
#' clamped to valid bounds.
#'
#' @param x A \code{MultiPhysioExperiment} object.
#' @param i Numeric vector of length 2 specifying time range \code{[tmin, tmax]}
#'   in seconds, or missing.
#' @param j Character vector of modality names, or missing.
#' @param ... Additional arguments (not used).
#' @param drop Logical (not used).
#' @return A \code{MultiPhysioExperiment} object with the selected
#'   modalities and/or time windows.
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(2500 * 4), nrow = 2500, ncol = 4)
#' emg_data <- matrix(rnorm(5000 * 2), nrow = 5000, ncol = 2)
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   samplingRate = 500
#' )
#' pe_emg <- PhysioExperiment(
#'   assays = list(raw = emg_data),
#'   samplingRate = 1000
#' )
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)
#'
#' # Subset by modality
#' mpe_eeg <- mpe[, "EEG"]
#' modalities(mpe_eeg)  # "EEG"
#'
#' # Subset by time range (seconds)
#' mpe_win <- mpe[c(1.0, 3.0), ]
#'
#' # Combined
#' mpe_sub <- mpe[c(1.0, 3.0), "EEG"]
setMethod("[", c("MultiPhysioExperiment", "ANY", "ANY"),
  function(x, i, j, ..., drop = FALSE) {
    exps <- experiments(x)
    mod_names <- names(exps)

    # ---- Modality subsetting (j) ----
    if (!missing(j)) {
      if (is.character(j)) {
        invalid <- setdiff(j, mod_names)
        if (length(invalid) > 0) {
          stop(
            "Unknown modality name(s): ",
            paste(invalid, collapse = ", "),
            call. = FALSE
          )
        }
        exps <- exps[j]
      } else {
        stop("j must be a character vector of modality names", call. = FALSE)
      }
    }

    # ---- Time-range subsetting (i) ----
    if (!missing(i)) {
      if (!is.numeric(i) || length(i) != 2) {
        stop(
          "i must be a numeric vector of length 2 specifying ",
          "time range [tmin, tmax] in seconds",
          call. = FALSE
        )
      }

      tmin <- i[1]
      tmax <- i[2]

      if (tmin >= tmax) {
        stop("tmin must be less than tmax", call. = FALSE)
      }

      exps <- lapply(exps, function(pe) {
        sr <- samplingRate(pe)
        n_samples <- length(pe)
        start_idx <- max(1L, as.integer(floor(tmin * sr)) + 1L)
        end_idx <- min(n_samples, as.integer(floor(tmax * sr)) + 1L)

        if (start_idx > n_samples || end_idx < 1L) {
          stop(
            "Time range [", tmin, ", ", tmax, "] is outside the data range ",
            "for modality with ", n_samples, " samples at ", sr, " Hz",
            call. = FALSE
          )
        }

        pe[start_idx:end_idx, ]
      })
    }

    # Rebuild the alignment for selected modalities
    new_alignment <- .build_default_alignment(exps)

    methods::new(
      "MultiPhysioExperiment",
      experiments = exps,
      alignment = new_alignment,
      sampleMap = S4Vectors::DataFrame(),
      couplingResults = list()
    )
  }
)
