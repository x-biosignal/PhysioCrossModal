# Directed coupling analysis: Granger causality
#
# Implements time-domain and spectral Granger causality for cross-modal
# physiological signal pairs. Fits bivariate autoregressive (AR) models
# and compares restricted vs unrestricted prediction variances.

#' Granger Causality Between Two Signals
#'
#' Computes directed coupling between two physiological signals using Granger
#' causality. The method fits bivariate autoregressive (AR) models to quantify
#' how much one signal helps predict the other beyond its own past. Both
#' time-domain and spectral (frequency-resolved) Granger causality are
#' supported.
#'
#' For the time-domain (parametric) method, Granger causality from x to y is
#' defined as:
#' \deqn{GC_{x \to y} = \log(\sigma^2_{restricted} / \sigma^2_{unrestricted})}
#' where the restricted model predicts y from its own past only, and the
#' unrestricted model predicts y from both its own past and x's past.
#'
#' For spectral Granger causality (when \code{freq_range} is specified), the AR
#' coefficients are transformed to the frequency domain via the transfer
#' function \eqn{H(f) = A(f)^{-1}}, and the spectral GC is computed as:
#' \deqn{GC_{x \to y}(f) = \log(S_{yy}^{restricted}(f) / S_{yy}^{unrestricted}(f))}
#' The returned values are then averaged over the specified frequency band.
#'
#' \strong{Note:} The spectral GC implementation uses separate restricted and
#' unrestricted univariate AR models, which is an approximation to the full
#' Geweke (1982) spectral decomposition based on a 2x2 bivariate VAR model.
#' Results may differ from toolboxes that implement the full Geweke decomposition
#' (e.g., MVGC, FieldTrip). For exact spectral GC, use the time-domain method
#' (\code{freq_range = NULL}) which correctly implements the Granger (1969)
#' log-ratio formulation.
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when \code{x} is a
#'   MultiPhysioExperiment.
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric
#'   vectors).
#' @param order Integer AR model order (default 5). Higher order captures
#'   longer-range temporal dependencies but requires more data.
#' @param freq_range Numeric vector of length 2 specifying the frequency band
#'   in Hz over which to compute spectral Granger causality. If \code{NULL}
#'   (default), returns time-domain GC values.
#' @param method Character; one of \code{"parametric"} (default) or
#'   \code{"nonparametric"}. Currently only parametric is implemented.
#' @param modality_x Character modality name in MultiPhysioExperiment for the
#'   x signal.
#' @param modality_y Character modality name in MultiPhysioExperiment for the
#'   y signal.
#' @param channels_x Integer channel index to extract from x (default 1).
#' @param channels_y Integer channel index to extract from y (default 1).
#'
#' @return A list with components:
#'   \describe{
#'     \item{gc_xy}{Numeric Granger causality from x to y (x drives y). A
#'       positive value indicates that x's past helps predict y beyond y's own
#'       past.}
#'     \item{gc_yx}{Numeric Granger causality from y to x (y drives x).}
#'     \item{net_gc}{Numeric net Granger causality (\code{gc_xy - gc_yx}).
#'       Positive values indicate that x drives y more than y drives x.}
#'     \item{order}{Integer AR model order used.}
#'   }
#'
#' @references
#' Granger, C. W. J. (1969). Investigating causal relations by econometric
#' models and cross-spectral methods. Econometrica, 37(3), 424-438.
#'
#' Geweke, J. (1982). Measurement of linear dependence and feedback between
#' multiple time series. Journal of the American Statistical Association,
#' 77(378), 304-313.
#'
#' @seealso [coherence()], [phaseLockingValue()]
#' @export
#' @examples
#' # Create directed signals: x drives y with a 5-sample lag
#' set.seed(42)
#' n <- 5000
#' x <- rnorm(n)
#' y <- numeric(n)
#' for (i in 6:n) y[i] <- 0.6 * x[i - 5] + 0.4 * rnorm(1)
#'
#' result <- grangerCausality(x, y, sr = 500, order = 10)
#' result$gc_xy   # Should be positive (x drives y)
#' result$net_gc  # Should be positive
grangerCausality <- function(x, y = NULL, sr = NULL, order = 5L,
                             freq_range = NULL,
                             method = c("parametric", "nonparametric"),
                             modality_x = NULL, modality_y = NULL,
                             channels_x = 1L, channels_y = 1L) {

  method <- match.arg(method)

  # Validate order
  order <- as.integer(order)
  if (length(order) != 1L || order < 1L) {
    stop("'order' must be a positive integer")
  }

  # Validate freq_range
  if (!is.null(freq_range)) {
    if (!is.numeric(freq_range) || length(freq_range) != 2 ||
        freq_range[1] >= freq_range[2]) {
      stop("'freq_range' must be a numeric vector of length 2 with ",
           "freq_range[1] < freq_range[2]")
    }
  }

  # Extract signal pair using shared utility
  pair <- .extract_signal_pair(
    x = x, y = y, sr = sr,
    modality_x = modality_x, modality_y = modality_y,
    channels_x = channels_x, channels_y = channels_y
  )
  sig_x <- pair$x
  sig_y <- pair$y
  sr_val <- pair$sr

  n <- length(sig_x)
  if (length(sig_y) != n) {
    stop("Signals x and y must have the same length after extraction")
  }
  if (n <= 2 * order + 1) {
    stop(sprintf(
      "Signal too short (%d samples) for AR model order %d. Need at least %d samples.",
      n, order, 2 * order + 2
    ))
  }

  # Demean signals
  sig_x <- sig_x - mean(sig_x)
  sig_y <- sig_y - mean(sig_y)

  # Compute GC in both directions
  if (is.null(freq_range)) {
    # Time-domain Granger causality
    gc_xy <- .gc_time_domain(sig_x, sig_y, order)
    gc_yx <- .gc_time_domain(sig_y, sig_x, order)
  } else {
    # Spectral Granger causality averaged over freq_range
    gc_xy <- .gc_spectral(sig_x, sig_y, order, sr_val, freq_range)
    gc_yx <- .gc_spectral(sig_y, sig_x, order, sr_val, freq_range)
  }

  list(
    gc_xy  = gc_xy,
    gc_yx  = gc_yx,
    net_gc = gc_xy - gc_yx,
    order  = order
  )
}


# ---- Internal Helpers --------------------------------------------------------

#' Compute time-domain Granger causality: does `driver` help predict `target`?
#'
#' Fits two AR models via OLS on lagged matrices:
#' - Restricted: target(t) ~ target(t-1) ... target(t-p)
#' - Unrestricted: target(t) ~ target(t-1) ... target(t-p), driver(t-1) ... driver(t-p)
#'
#' GC = log(var_restricted / var_unrestricted), clamped to >= 0.
#'
#' @param driver Numeric vector (demeaned).
#' @param target Numeric vector (demeaned, same length as driver).
#' @param order Integer AR model order.
#' @return Numeric scalar GC value (>= 0).
#' @noRd
.gc_time_domain <- function(driver, target, order) {
  n <- length(target)
  n_obs <- n - order  # number of usable observations

  # Build the response vector: target values from (order+1) to n
  y_vec <- target[(order + 1L):n]

  # Build lagged matrix for restricted model (target lags only)
  X_restricted <- matrix(0, nrow = n_obs, ncol = order)
  for (k in seq_len(order)) {
    X_restricted[, k] <- target[(order + 1L - k):(n - k)]
  }

  # Build lagged matrix for unrestricted model (target + driver lags)
  X_unrestricted <- matrix(0, nrow = n_obs, ncol = 2L * order)
  for (k in seq_len(order)) {
    X_unrestricted[, k] <- target[(order + 1L - k):(n - k)]
    X_unrestricted[, order + k] <- driver[(order + 1L - k):(n - k)]
  }

  # Fit via manual OLS: beta = solve(X'X) %*% X'y
  # Compute residual variance for each model
  var_restricted <- .ols_residual_var(X_restricted, y_vec)
  var_unrestricted <- .ols_residual_var(X_unrestricted, y_vec)

  # GC = log(var_restricted / var_unrestricted), clamped to >= 0
  if (var_unrestricted < .Machine$double.eps ||
      var_restricted < .Machine$double.eps) {
    return(0)
  }
  max(0, log(var_restricted / var_unrestricted))
}


#' OLS residual variance with ridge regularization
#'
#' Solves y = X * beta via normal equations with small ridge penalty for
#' numerical stability, and returns the variance of the residuals.
#'
#' @param X Numeric matrix (n_obs x p).
#' @param y Numeric vector (n_obs).
#' @return Numeric scalar residual variance.
#' @noRd
.ols_residual_var <- function(X, y) {
  XtX <- crossprod(X)          # X'X  (p x p)
  Xty <- crossprod(X, y)       # X'y  (p x 1)

  # Ridge regularization for numerical stability
  reg_lambda <- 1e-8 * max(abs(diag(XtX)))
  diag(XtX) <- diag(XtX) + reg_lambda

  beta <- tryCatch(
    solve(XtX, Xty),
    error = function(e) rep(0, ncol(X))
  )

  residuals <- as.numeric(y - X %*% beta)
  mean(residuals^2)
}


#' Compute spectral Granger causality: does `driver` help predict `target`?
#'
#' Fits AR models via OLS, then transforms to frequency domain via transfer
#' functions. GC(f) = log(S_restricted(f) / S_unrestricted(f)) averaged over
#' the specified frequency range.
#'
#' @param driver Numeric vector (demeaned).
#' @param target Numeric vector (demeaned, same length as driver).
#' @param order Integer AR model order.
#' @param sr Numeric sampling rate in Hz.
#' @param freq_range Numeric vector of length 2 (low, high) in Hz.
#' @return Numeric scalar spectral GC value averaged over freq_range (>= 0).
#' @noRd
.gc_spectral <- function(driver, target, order, sr, freq_range) {
  n <- length(target)
  n_obs <- n - order

  # Build response
  y_vec <- target[(order + 1L):n]

  # ---- Restricted model (target lags only) ----
  X_restricted <- matrix(0, nrow = n_obs, ncol = order)
  for (k in seq_len(order)) {
    X_restricted[, k] <- target[(order + 1L - k):(n - k)]
  }
  ar_restricted <- .ols_coefficients(X_restricted, y_vec)
  var_restricted <- .ols_residual_var(X_restricted, y_vec)

  # ---- Unrestricted model (target + driver lags) ----
  X_unrestricted <- matrix(0, nrow = n_obs, ncol = 2L * order)
  for (k in seq_len(order)) {
    X_unrestricted[, k] <- target[(order + 1L - k):(n - k)]
    X_unrestricted[, order + k] <- driver[(order + 1L - k):(n - k)]
  }
  ar_unrestricted <- .ols_coefficients(X_unrestricted, y_vec)
  var_unrestricted <- .ols_residual_var(X_unrestricted, y_vec)

  # Extract coefficients for target-lags and driver-lags from unrestricted model
  a_target <- ar_unrestricted[seq_len(order)]
  a_driver <- ar_unrestricted[order + seq_len(order)]

  # ---- Compute spectral GC over frequency grid ----
  n_freq <- 128L
  freqs <- seq(0, sr / 2, length.out = n_freq)

  gc_spectrum <- numeric(n_freq)

  for (fi in seq_len(n_freq)) {
    f <- freqs[fi]

    # Restricted transfer function: H_r(f) = 1 - sum(a_k * exp(-2*pi*i*f*k/sr))
    H_restricted <- 1 + 0i
    for (k in seq_len(order)) {
      H_restricted <- H_restricted - ar_restricted[k] * exp(-2i * pi * f * k / sr)
    }

    # Unrestricted transfer function for target-only part
    H_target <- 1 + 0i
    for (k in seq_len(order)) {
      H_target <- H_target - a_target[k] * exp(-2i * pi * f * k / sr)
    }

    # Spectral densities from restricted and unrestricted models
    S_restricted <- var_restricted / Mod(H_restricted)^2
    S_unrestricted <- var_unrestricted / Mod(H_target)^2

    if (S_unrestricted > .Machine$double.eps &&
        S_restricted > .Machine$double.eps) {
      gc_spectrum[fi] <- max(0, log(S_restricted / S_unrestricted))
    }
  }

  # Average over requested frequency range
  band_idx <- which(freqs >= freq_range[1] & freqs <= freq_range[2])
  if (length(band_idx) == 0L) {
    return(0)
  }
  max(0, mean(gc_spectrum[band_idx]))
}


#' OLS coefficient estimation with ridge regularization
#'
#' @param X Numeric matrix (n_obs x p).
#' @param y Numeric vector (n_obs).
#' @return Numeric vector of estimated coefficients (length p).
#' @noRd
.ols_coefficients <- function(X, y) {
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)

  reg_lambda <- 1e-8 * max(abs(diag(XtX)))
  diag(XtX) <- diag(XtX) + reg_lambda

  beta <- tryCatch(
    solve(XtX, Xty),
    error = function(e) rep(0, ncol(X))
  )
  as.numeric(beta)
}
