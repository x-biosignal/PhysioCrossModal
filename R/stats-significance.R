# Statistical significance testing for coupling analysis
#
# Provides surrogate testing and bootstrap confidence intervals for
# all coupling methods in PhysioCrossModal.

# ---- Internal surrogate generators ------------------------------------------

#' Phase-randomization surrogate
#'
#' Randomises the FFT phases of a signal while preserving its amplitude
#' spectrum (and thus its autocorrelation structure). Hermitian symmetry is
#' enforced so the result is real-valued.
#'
#' @param x Numeric vector.
#' @return Numeric vector of the same length with randomised phases.
#' @noRd
.phase_randomize_surrogate <- function(x) {
  n <- length(x)
  fft_x <- stats::fft(x)
  amp <- Mod(fft_x)

  # Generate random phases with Hermitian symmetry
  if (n %% 2 == 0) {
    half <- n / 2
    # DC (index 1) and Nyquist (index half+1) must have zero phase
    rand_phase <- stats::runif(half - 1L, -pi, pi)
    phases <- c(0, rand_phase, 0, -rev(rand_phase))
  } else {
    half <- (n - 1) / 2
    rand_phase <- stats::runif(half, -pi, pi)
    phases <- c(0, rand_phase, -rev(rand_phase))
  }

  fft_surr <- amp * exp(complex(imaginary = phases))
  Re(stats::fft(fft_surr, inverse = TRUE)) / n
}

#' Time-shift surrogate
#'
#' Circularly shifts a signal by a random offset, destroying the temporal
#' relationship with any other signal while preserving autocorrelation
#' structure exactly.
#'
#' @param x Numeric vector.
#' @return Numeric vector of the same length, circularly shifted.
#' @noRd
.timeshift_surrogate <- function(x) {
  n <- length(x)
  shift <- sample.int(n - 1L, 1L)
  c(x[(shift + 1L):n], x[seq_len(shift)])
}

#' Moving block bootstrap indices
#'
#' Generates indices for a moving-block bootstrap sample. Blocks of
#' consecutive indices are randomly drawn (with replacement) and
#' concatenated to form a bootstrap sample of length \code{n}.
#'
#' @param n Integer sample size.
#' @param block_len Integer block length. Defaults to \code{ceiling(sqrt(n))}.
#' @return Integer vector of length \code{n} with bootstrap indices.
#' @noRd
.block_bootstrap_indices <- function(n, block_len = NULL) {
  if (is.null(block_len)) block_len <- ceiling(sqrt(n))
  block_len <- max(1L, as.integer(block_len))

  n_blocks <- ceiling(n / block_len)
  max_start <- n - block_len + 1L

  starts <- sample.int(max_start, n_blocks, replace = TRUE)
  indices <- unlist(lapply(starts, function(s) s:(s + block_len - 1L)))
  indices[seq_len(n)]
}

# ---- Parallel computation helper ---------------------------------------------

#' Apply a function in parallel or sequentially
#'
#' Uses \code{parallel::mclapply()} on Unix systems and
#' \code{parallel::parLapply()} on Windows. Falls back to sequential
#' \code{lapply()} on error.
#'
#' @param X A vector to iterate over.
#' @param FUN Function to apply to each element.
#' @param cores Integer number of cores.
#' @return A list of results.
#' @noRd
.parallel_apply <- function(X, FUN, cores) {
  tryCatch({
    if (.Platform$OS.type == "unix") {
      parallel::mclapply(X, FUN, mc.cores = cores)
    } else {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      parallel::parLapply(cl, X, FUN)
    }
  }, error = function(e) {
    lapply(X, FUN)
  })
}

# ---- Exported functions ------------------------------------------------------

#' Surrogate-based significance test for coupling
#'
#' Tests whether the observed coupling between two signals is statistically
#' significant by comparing it against a null distribution generated from
#' surrogate signals. Surrogates are created by randomising the phases (or
#' circularly shifting) one signal, thereby destroying the cross-signal
#' relationship while preserving the autocorrelation structure.
#'
#' The p-value is computed as \code{(sum(surr >= obs) + 1) / (n_surr + 1)},
#' following the conservative correction of Phipson & Smyth (2010).
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when \code{x} is an
#'   MPE.
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param method Character coupling method, one of \code{"coherence"},
#'   \code{"plv"}, \code{"pli"}, \code{"wpli"}, \code{"granger"},
#'   \code{"crosscorrelation"}, \code{"wavelet_coherence"},
#'   \code{"wavelet_plv"}.
#' @param n_surrogates Integer number of surrogates (default 199).
#' @param surrogate_type Character surrogate generation method:
#'   \code{"phase"} (default) or \code{"timeshift"}.
#' @param modality_x,modality_y Character modality names for MPE input.
#' @param channels_x,channels_y Integer channel indices (default 1).
#' @param cores Integer number of parallel cores to use (default \code{1L}).
#'   When \code{cores > 1}, surrogate computations are distributed across
#'   cores using \code{parallel::mclapply()} (Unix) or
#'   \code{parallel::parLapply()} (Windows). Falls back to sequential
#'   computation on error.
#' @param ... Additional arguments passed to the coupling function
#'   (e.g. \code{freq_band}, \code{nperseg}, \code{order}).
#'
#' @return A list with components:
#'   \describe{
#'     \item{observed}{The full coupling result from the original signals.}
#'     \item{statistic}{Numeric scalar: the extracted coupling statistic.}
#'     \item{p_value}{Numeric scalar: surrogate p-value.}
#'     \item{surrogate_distribution}{Numeric vector of surrogate statistics.}
#'     \item{threshold_95}{Numeric scalar: 95th percentile of surrogates.}
#'   }
#'
#' @references
#' Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., & Farmer, J. D.
#' (1992). Testing for nonlinearity in time series: the method of surrogate
#' data. \emph{Physica D: Nonlinear Phenomena}, 58(1--4), 77--94.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never
#' be zero: calculating exact P-values when permutations are randomly drawn.
#' \emph{Statistical Applications in Genetics and Molecular Biology}, 9(1),
#' Article 39.
#'
#' @seealso [bootstrapCI()], [surrogateMatrixTest()],
#'   [couplingAnalysis()]
#' @export
#' @examples
#' sr <- 500
#' t <- seq(0, 2, length.out = sr * 2)
#' x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
#' y <- 0.8 * sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
#' result <- surrogateTest(x, y, sr = sr, method = "coherence",
#'                         n_surrogates = 19, nperseg = 128L)
#' result$p_value
surrogateTest <- function(x, y = NULL, sr = NULL, method,
                          n_surrogates = 199L,
                          surrogate_type = c("phase", "timeshift"),
                          modality_x = NULL, modality_y = NULL,
                          channels_x = 1L, channels_y = 1L,
                          cores = 1L,
                          ...) {

  surrogate_type <- match.arg(surrogate_type)
  n_surrogates <- as.integer(n_surrogates)
  if (n_surrogates < 1L) stop("n_surrogates must be >= 1", call. = FALSE)

  # Extract signal pair once
  pair <- .extract_signal_pair(x, y, sr = sr,
                               modality_x = modality_x,
                               modality_y = modality_y,
                               channels_x = channels_x,
                               channels_y = channels_y)
  sig_x <- pair$x
  sig_y <- pair$y
  sr_val <- pair$sr

  # Compute observed coupling
  observed <- .dispatch_coupling(sig_x, sig_y, sr_val, method, ...)
  obs_stat <- .extract_coupling_stat(observed, method)

  # Generate surrogate distribution
  surr_gen <- switch(surrogate_type,
    phase = .phase_randomize_surrogate,
    timeshift = .timeshift_surrogate
  )

  cores <- as.integer(cores)
  if (cores > 1L) {
    surr_stats <- unlist(.parallel_apply(seq_len(n_surrogates), function(i) {
      y_surr <- surr_gen(sig_y)
      surr_result <- .dispatch_coupling(sig_x, y_surr, sr_val, method, ...)
      .extract_coupling_stat(surr_result, method)
    }, cores = cores))
  } else {
    surr_stats <- numeric(n_surrogates)
    for (i in seq_len(n_surrogates)) {
      y_surr <- surr_gen(sig_y)
      surr_result <- .dispatch_coupling(sig_x, y_surr, sr_val, method, ...)
      surr_stats[i] <- .extract_coupling_stat(surr_result, method)
    }
  }

  # Conservative p-value (Phipson & Smyth 2010)
  p_value <- (sum(surr_stats >= obs_stat) + 1L) / (n_surrogates + 1L)

  list(
    observed = observed,
    statistic = obs_stat,
    p_value = p_value,
    surrogate_distribution = surr_stats,
    threshold_95 = as.numeric(stats::quantile(surr_stats, 0.95))
  )
}


#' Bootstrap confidence interval for coupling
#'
#' Computes a bootstrap confidence interval for the coupling statistic
#' between two signals using the moving-block bootstrap (to preserve
#' temporal autocorrelation).
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when \code{x} is an
#'   MPE.
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param method Character coupling method (same options as
#'   \code{\link{surrogateTest}}).
#' @param n_boot Integer number of bootstrap replicates (default 199).
#' @param ci Numeric confidence level (default 0.95).
#' @param block_len Integer block length for block bootstrap, or NULL for
#'   automatic (\code{ceiling(sqrt(n))}).
#' @param modality_x,modality_y Character modality names for MPE input.
#' @param channels_x,channels_y Integer channel indices (default 1).
#' @param cores Integer number of parallel cores to use (default \code{1L}).
#'   When \code{cores > 1}, bootstrap replicates are computed in parallel.
#' @param ... Additional arguments passed to the coupling function.
#'
#' @return A list with components:
#'   \describe{
#'     \item{observed}{The full coupling result from the original signals.}
#'     \item{statistic}{Numeric scalar: the extracted coupling statistic.}
#'     \item{ci_lower}{Numeric scalar: lower CI bound.}
#'     \item{ci_upper}{Numeric scalar: upper CI bound.}
#'     \item{ci_level}{Numeric scalar: confidence level used.}
#'     \item{boot_distribution}{Numeric vector of bootstrap statistics.}
#'   }
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An Introduction to the
#' Bootstrap}. Chapman & Hall/CRC.
#'
#' @seealso [surrogateTest()], [couplingAnalysis()]
#' @export
#' @examples
#' sr <- 500
#' t <- seq(0, 2, length.out = sr * 2)
#' x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
#' y <- 0.8 * sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
#' result <- bootstrapCI(x, y, sr = sr, method = "coherence",
#'                       n_boot = 19, nperseg = 128L)
#' c(result$ci_lower, result$ci_upper)
bootstrapCI <- function(x, y = NULL, sr = NULL, method,
                        n_boot = 199L, ci = 0.95,
                        block_len = NULL,
                        modality_x = NULL, modality_y = NULL,
                        channels_x = 1L, channels_y = 1L,
                        cores = 1L,
                        ...) {

  n_boot <- as.integer(n_boot)
  if (n_boot < 1L) stop("n_boot must be >= 1", call. = FALSE)
  if (ci <= 0 || ci >= 1) stop("ci must be between 0 and 1", call. = FALSE)

  # Extract signal pair once
  pair <- .extract_signal_pair(x, y, sr = sr,
                               modality_x = modality_x,
                               modality_y = modality_y,
                               channels_x = channels_x,
                               channels_y = channels_y)
  sig_x <- pair$x
  sig_y <- pair$y
  sr_val <- pair$sr

  n <- min(length(sig_x), length(sig_y))
  sig_x <- sig_x[seq_len(n)]
  sig_y <- sig_y[seq_len(n)]

  # Compute observed coupling
  observed <- .dispatch_coupling(sig_x, sig_y, sr_val, method, ...)
  obs_stat <- .extract_coupling_stat(observed, method)

  # Block bootstrap
  cores <- as.integer(cores)
  if (cores > 1L) {
    boot_stats <- unlist(.parallel_apply(seq_len(n_boot), function(i) {
      idx <- .block_bootstrap_indices(n, block_len)
      boot_result <- .dispatch_coupling(sig_x[idx], sig_y[idx], sr_val,
                                        method, ...)
      .extract_coupling_stat(boot_result, method)
    }, cores = cores))
  } else {
    boot_stats <- numeric(n_boot)
    for (i in seq_len(n_boot)) {
      idx <- .block_bootstrap_indices(n, block_len)
      boot_result <- .dispatch_coupling(sig_x[idx], sig_y[idx], sr_val,
                                        method, ...)
      boot_stats[i] <- .extract_coupling_stat(boot_result, method)
    }
  }

  alpha <- 1 - ci
  ci_lower <- as.numeric(stats::quantile(boot_stats, alpha / 2))
  ci_upper <- as.numeric(stats::quantile(boot_stats, 1 - alpha / 2))

  list(
    observed = observed,
    statistic = obs_stat,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    ci_level = ci,
    boot_distribution = boot_stats
  )
}

# ---- Surrogate matrix test ---------------------------------------------------

#' Surrogate-based significance test for coupling matrices
#'
#' Tests each element of a coupling matrix for significance using
#' surrogate testing, with correction for multiple comparisons
#' (FDR or Bonferroni).
#'
#' @param x PhysioExperiment or MultiPhysioExperiment.
#' @param y PhysioExperiment or NULL (when \code{x} is an MPE).
#' @param method Character coupling method (same options as
#'   \code{\link{couplingAnalysis}}).
#' @param n_surrogates Integer number of surrogates (default 199).
#' @param surrogate_type Character surrogate generation method:
#'   \code{"phase"} (default) or \code{"timeshift"}.
#' @param correction Character correction method: \code{"fdr"} (default),
#'   \code{"bonferroni"}, or \code{"none"}.
#' @param alpha Numeric significance level (default 0.05).
#' @param channels_x,channels_y Integer vectors of channel indices, or
#'   NULL for all channels.
#' @param modality_x,modality_y Character modality names for MPE input.
#' @param cores Integer number of parallel cores (default 1).
#' @param ... Additional arguments passed to the coupling function.
#'
#' @return A list with components:
#'   \describe{
#'     \item{matrix}{Coupling matrix (observed values).}
#'     \item{p_values}{Matrix of raw p-values (same dimensions).}
#'     \item{p_adjusted}{Matrix of corrected p-values.}
#'     \item{significant}{Logical matrix indicating significance.}
#'     \item{correction}{Character correction method used.}
#'     \item{alpha}{Numeric significance level.}
#'   }
#'
#' @references
#' Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., & Farmer, J. D.
#' (1992). Testing for nonlinearity in time series: the method of surrogate
#' data. \emph{Physica D: Nonlinear Phenomena}, 58(1--4), 77--94.
#'
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing. \emph{Journal
#' of the Royal Statistical Society: Series B (Methodological)}, 57(1),
#' 289--300.
#'
#' @seealso [surrogateTest()], [couplingMatrix()], [coherenceMatrix()]
#' @export
surrogateMatrixTest <- function(x, y = NULL, method,
                                 n_surrogates = 199L,
                                 surrogate_type = c("phase", "timeshift"),
                                 correction = c("fdr", "bonferroni", "none"),
                                 alpha = 0.05,
                                 channels_x = NULL, channels_y = NULL,
                                 modality_x = NULL, modality_y = NULL,
                                 cores = 1L, ...) {

  surrogate_type <- match.arg(surrogate_type)
  correction <- match.arg(correction)
  n_surrogates <- as.integer(n_surrogates)

  # Resolve PE pair
  pe_pair <- .resolve_pe_pair(x, y, modality_x, modality_y)
  pe_x <- pe_pair$pe_x
  pe_y <- pe_pair$pe_y

  # Resolve channel indices
  n_ch_x <- ncol(SummarizedExperiment::assay(pe_x, 1L))
  n_ch_y <- ncol(SummarizedExperiment::assay(pe_y, 1L))
  if (is.null(channels_x)) channels_x <- seq_len(n_ch_x)
  if (is.null(channels_y)) channels_y <- seq_len(n_ch_y)

  # Channel names
  cd_x <- SummarizedExperiment::colData(pe_x)
  cd_y <- SummarizedExperiment::colData(pe_y)
  ch_names_x <- if ("label" %in% colnames(cd_x)) {
    as.character(cd_x$label[channels_x])
  } else {
    paste0("x_ch", channels_x)
  }
  ch_names_y <- if ("label" %in% colnames(cd_y)) {
    as.character(cd_y$label[channels_y])
  } else {
    paste0("y_ch", channels_y)
  }

  # Compute observed coupling matrix and p-values
  obs_mat <- matrix(NA_real_, nrow = length(channels_x),
                    ncol = length(channels_y),
                    dimnames = list(ch_names_x, ch_names_y))
  p_mat <- matrix(NA_real_, nrow = length(channels_x),
                  ncol = length(channels_y),
                  dimnames = list(ch_names_x, ch_names_y))

  for (i in seq_along(channels_x)) {
    for (j in seq_along(channels_y)) {
      result <- surrogateTest(
        pe_x, pe_y,
        method = method,
        n_surrogates = n_surrogates,
        surrogate_type = surrogate_type,
        channels_x = channels_x[i],
        channels_y = channels_y[j],
        cores = cores,
        ...
      )
      obs_mat[i, j] <- result$statistic
      p_mat[i, j] <- result$p_value
    }
  }

  # Apply multiple comparison correction
  p_vec <- as.vector(p_mat)
  p_adj_vec <- switch(correction,
    fdr = stats::p.adjust(p_vec, method = "fdr"),
    bonferroni = stats::p.adjust(p_vec, method = "bonferroni"),
    none = p_vec
  )
  p_adj <- matrix(p_adj_vec, nrow = nrow(p_mat), ncol = ncol(p_mat),
                  dimnames = dimnames(p_mat))

  sig_mat <- p_adj < alpha

  list(
    matrix = obs_mat,
    p_values = p_mat,
    p_adjusted = p_adj,
    significant = sig_mat,
    correction = correction,
    alpha = alpha
  )
}

# ---- Internal dispatch -------------------------------------------------------

#' Dispatch coupling computation for raw numeric vectors
#'
#' @param sig_x Numeric vector.
#' @param sig_y Numeric vector.
#' @param sr Numeric sampling rate.
#' @param method Character coupling method.
#' @param ... Additional arguments.
#' @return Coupling result list.
#' @noRd
.dispatch_coupling <- function(sig_x, sig_y, sr, method, ...) {
  switch(method,
    coherence = coherence(sig_x, sig_y, sr = sr, ...),
    plv = phaseLockingValue(sig_x, sig_y, sr = sr, ...),
    pli = phaseLagIndex(sig_x, sig_y, sr = sr, ...),
    wpli = weightedPLI(sig_x, sig_y, sr = sr, ...),
    granger = grangerCausality(sig_x, sig_y, sr = sr, ...),
    crosscorrelation = crossCorrelation(sig_x, sig_y, sr = sr, ...),
    wavelet_coherence = waveletCoherence(sig_x, sig_y, sr = sr, ...),
    wavelet_plv = waveletPLV(sig_x, sig_y, sr = sr, ...),
    multitaper_coherence = multitaperCoherence(sig_x, sig_y, sr = sr, ...),
    stop("Unknown coupling method: '", method, "'", call. = FALSE)
  )
}
