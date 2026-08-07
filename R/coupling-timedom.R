# Time-domain coupling functions for cross-modal analysis
#
# Functions for computing cross-correlation between signals from different
# physiological modalities, including both full and sliding-window variants.

#' Cross-correlation between two signals
#'
#' Computes the cross-correlation between two signals at various lags,
#' measuring their time-domain coupling. Cross-correlation quantifies how
#' similar one signal is to a time-shifted version of another. A positive
#' peak lag indicates that `y` leads `x` (i.e., `y` must be shifted forward
#' to align with `x`), while a negative peak lag indicates `x` leads `y`.
#'
#' The function accepts numeric vectors, PhysioExperiment objects, or a
#' MultiPhysioExperiment with named modalities, using
#' \code{.extract_signal_pair()} internally for flexible input handling.
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when \code{x} is
#'   a MultiPhysioExperiment.
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param max_lag Integer maximum lag in samples. If NULL, defaults to
#'   \code{floor(length(x) / 4)}.
#' @param normalize Logical; if TRUE (default), compute normalized
#'   cross-correlation (Pearson-like, values in \[-1, 1\]).
#' @param modality_x Character modality name in MPE for the x signal.
#' @param modality_y Character modality name in MPE for the y signal.
#' @param channels_x Integer which channel to extract from x (default 1).
#' @param channels_y Integer which channel to extract from y (default 1).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{correlation}{Numeric vector of cross-correlation values at each
#'       lag.}
#'     \item{lags}{Integer vector of lag values in samples.}
#'     \item{lag_seconds}{Numeric vector of lag values in seconds.}
#'     \item{peak_lag}{Integer lag (in samples) at which the absolute
#'       correlation is maximised.}
#'     \item{peak_lag_seconds}{Numeric peak lag converted to seconds.}
#'     \item{peak_correlation}{Numeric cross-correlation value at the peak
#'       lag.}
#'   }
#'
#' @references
#' Chatfield, C. (2004). \emph{The Analysis of Time Series: An
#' Introduction} (6th ed.). Chapman & Hall/CRC.
#'
#' @seealso [slidingCrossCorrelation()], [coherence()],
#'   [couplingAnalysis()]
#' @export
#' @examples
#' # Simple example: y is a delayed copy of x
#' sr <- 500
#' set.seed(1)
#' x <- rnorm(1000)
#' y <- c(rep(0, 10), x[1:990])
#' result <- crossCorrelation(x, y, sr = sr)
#' result$peak_lag       # should be 10
#' result$peak_lag_seconds
crossCorrelation <- function(x, y = NULL, sr = NULL, max_lag = NULL,
                             normalize = TRUE,
                             modality_x = NULL, modality_y = NULL,
                             channels_x = 1L, channels_y = 1L) {

  # Extract signal pair from various input types
  pair <- .extract_signal_pair(
    x, y, sr = sr,
    modality_x = modality_x, modality_y = modality_y,
    channels_x = channels_x, channels_y = channels_y
  )

  sig_x <- pair$x
  sig_y <- pair$y
  sr_val <- pair$sr

  # Ensure equal length (trim to shorter if needed)
  min_len <- min(length(sig_x), length(sig_y))
  sig_x <- sig_x[seq_len(min_len)]
  sig_y <- sig_y[seq_len(min_len)]

  n <- length(sig_x)

  if (n < 2L) {
    stop("Signals must have at least 2 samples for cross-correlation",
         call. = FALSE)
  }

  # Default max_lag to length / 4

if (is.null(max_lag)) {
    max_lag <- floor(n / 4L)
  }
  max_lag <- as.integer(max_lag)

  if (max_lag < 0L) {
    stop("max_lag must be non-negative", call. = FALSE)
  }
  if (max_lag >= n) {
    max_lag <- n - 1L
  }

  # Compute cross-correlation using stats::ccf
  # stats::ccf returns normalized (Pearson-like) values by default
  # type = "correlation" normalizes, type = "covariance" does not
  ccf_type <- if (normalize) "correlation" else "covariance"
  ccf_result <- stats::ccf(sig_x, sig_y, lag.max = max_lag, plot = FALSE,
                            type = ccf_type)

  lags <- as.integer(ccf_result$lag[, 1, 1])
  correlation <- as.numeric(ccf_result$acf[, 1, 1])
  lag_seconds <- lags / sr_val

  # Find peak: the lag with maximum absolute correlation
  peak_idx <- which.max(abs(correlation))
  peak_lag <- lags[peak_idx]
  peak_lag_seconds <- peak_lag / sr_val
  peak_correlation <- correlation[peak_idx]

  list(
    correlation = correlation,
    lags = lags,
    lag_seconds = lag_seconds,
    peak_lag = peak_lag,
    peak_lag_seconds = peak_lag_seconds,
    peak_correlation = peak_correlation
  )
}


#' Sliding-window cross-correlation
#'
#' Computes cross-correlation in sliding (overlapping) windows to track
#' how time-domain coupling varies over time. For each window position,
#' \code{\link{crossCorrelation}} is called and the results are assembled
#' into a matrix.
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when \code{x} is
#'   a MultiPhysioExperiment.
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param window_sec Numeric window length in seconds (default 1).
#' @param step_sec Numeric step size in seconds (default 0.5).
#' @param max_lag Integer maximum lag in samples for each window. If NULL,
#'   defaults to \code{floor(window_samples / 4)}.
#' @param modality_x Character modality name in MPE for the x signal.
#' @param modality_y Character modality name in MPE for the y signal.
#' @param channels_x Integer which channel to extract from x (default 1).
#' @param channels_y Integer which channel to extract from y (default 1).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{correlations}{Numeric matrix of dimensions (n_windows x n_lags)
#'       containing cross-correlation values.}
#'     \item{times}{Numeric vector of window centre times in seconds.}
#'     \item{lags}{Integer vector of lag values in samples.}
#'     \item{peak_lags}{Numeric vector of peak lags (in samples) for each
#'       window.}
#'     \item{peak_correlations}{Numeric vector of peak correlation values
#'       for each window.}
#'   }
#'
#' @references
#' Chatfield, C. (2004). \emph{The Analysis of Time Series: An
#' Introduction} (6th ed.). Chapman & Hall/CRC.
#'
#' @seealso [crossCorrelation()], [plotCouplingTimecourse()],
#'   [couplingAnalysis()]
#' @export
#' @examples
#' sr <- 500
#' set.seed(1)
#' x <- rnorm(5000)
#' y <- c(rep(0, 10), x[1:4990])
#' result <- slidingCrossCorrelation(x, y, sr = sr,
#'                                    window_sec = 1, step_sec = 0.5)
#' dim(result$correlations)
#' result$peak_lags
slidingCrossCorrelation <- function(x, y = NULL, sr = NULL,
                                     window_sec = 1, step_sec = 0.5,
                                     max_lag = NULL,
                                     modality_x = NULL, modality_y = NULL,
                                     channels_x = 1L, channels_y = 1L) {

  # Extract signal pair from various input types
  pair <- .extract_signal_pair(
    x, y, sr = sr,
    modality_x = modality_x, modality_y = modality_y,
    channels_x = channels_x, channels_y = channels_y
  )

  sig_x <- pair$x
  sig_y <- pair$y
  sr_val <- pair$sr

  # Ensure equal length
  min_len <- min(length(sig_x), length(sig_y))
  sig_x <- sig_x[seq_len(min_len)]
  sig_y <- sig_y[seq_len(min_len)]

  n <- length(sig_x)
  duration <- n / sr_val

  # Convert window/step from seconds to samples
  window_samples <- as.integer(round(window_sec * sr_val))
  step_samples <- as.integer(round(step_sec * sr_val))

  if (window_samples < 2L) {
    stop("Window size must be at least 2 samples", call. = FALSE)
  }
  if (step_samples < 1L) {
    stop("Step size must be at least 1 sample", call. = FALSE)
  }
  if (window_samples > n) {
    stop(sprintf(
      "Window size (%d samples = %.2f s) exceeds signal length (%d samples = %.2f s)",
      window_samples, window_sec, n, duration
    ), call. = FALSE)
  }

  # Default max_lag for each window
  if (is.null(max_lag)) {
    max_lag <- floor(window_samples / 4L)
  }
  max_lag <- as.integer(max_lag)

  # Compute window positions
  n_windows <- floor((n - window_samples) / step_samples) + 1L
  if (n_windows < 1L) {
    stop("Signal too short for the specified window and step sizes",
         call. = FALSE)
  }

  # Pre-compute lag vector from first window to allocate output
  lags <- seq(-max_lag, max_lag)
  n_lags <- length(lags)

  # Allocate output
  correlations <- matrix(NA_real_, nrow = n_windows, ncol = n_lags)
  times <- numeric(n_windows)
  peak_lags <- numeric(n_windows)
  peak_correlations <- numeric(n_windows)

  for (w in seq_len(n_windows)) {
    start <- (w - 1L) * step_samples + 1L
    end <- start + window_samples - 1L

    win_x <- sig_x[start:end]
    win_y <- sig_y[start:end]

    # Compute cross-correlation for this window
    cc <- crossCorrelation(win_x, win_y, sr = sr_val,
                           max_lag = max_lag, normalize = TRUE)

    correlations[w, ] <- cc$correlation
    times[w] <- ((start - 1L) + (end - 1L)) / 2 / sr_val  # window centre
    peak_lags[w] <- cc$peak_lag
    peak_correlations[w] <- cc$peak_correlation
  }

  list(
    correlations = correlations,
    times = times,
    lags = lags,
    peak_lags = peak_lags,
    peak_correlations = peak_correlations
  )
}
