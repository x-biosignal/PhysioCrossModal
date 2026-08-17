# Directed coupling analysis: Granger causality
#
# Implements time-domain and spectral Granger causality for cross-modal
# physiological signal pairs. Fits bivariate autoregressive (AR) models
# and compares restricted vs unrestricted prediction variances.

#' Granger Causality Between Two Signals
#'
#' Computes directed coupling between two physiological signals using either
#' autoregressive prediction or nonparametric spectral factorization. Both
#' scalar and frequency-resolved Granger causality are supported.
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
#' The parametric frequency-band path uses separate restricted and unrestricted
#' univariate AR models and is an approximation to the full Geweke spectral
#' decomposition. Its established output is retained for compatibility. The
#' nonparametric path instead estimates a complete Hermitian 2x2
#' cross-spectral matrix with common windows/tapers, applies Wilson's
#' minimum-phase factorization, and evaluates Geweke directional spectra. A
#' scalar is the trapezoidal frequency mean over the selected band, with DC
#' and Nyquist receiving endpoint half weights.
#'
#' Cross-spectral eigenvalue adjustment is limited to the reported relative
#' `eigen_floor`; materially indefinite matrices fail. Wilson update
#' convergence and spectral reconstruction are checked separately. The
#' relative reconstruction limit is
#' `min(0.05, max(0.02, 100 * factor_tolerance))`, so it remains between 2%
#' and 5%; the observed error and limit are returned in the result.
#' Nonparametric estimation requires at least 32 paired samples and at least
#' two spectral realizations.
#'
#' Granger causality is a stationary linear predictive-direction measure. It
#' is not proof of physical, mechanistic, or interventional causality and does
#' not resolve unmeasured confounding, common drivers, volume conduction,
#' instantaneous mixing, or non-stationarity.
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when \code{x} is a
#'   MultiPhysioExperiment.
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric
#'   vectors).
#' @param order Integer AR model order (default 5). Higher order captures
#'   longer-range temporal dependencies but requires more data.
#' @param freq_range Numeric vector of length 2 specifying the frequency band
#'   in Hz over which to compute spectral Granger causality. With the
#'   parametric method, `NULL` returns time-domain GC. With the nonparametric
#'   method, `NULL` aggregates the full DC-to-Nyquist spectrum.
#' @param method Character estimator. `"parametric"` (default) uses the
#'   established MVAR ordinary-least-squares path. `"nonparametric"` estimates
#'   a complete cross-spectral matrix and applies Wilson factorization followed
#'   by Geweke's frequency-domain decomposition.
#' @param modality_x Character modality name in MultiPhysioExperiment for the
#'   x signal.
#' @param modality_y Character modality name in MultiPhysioExperiment for the
#'   y signal.
#' @param channels_x Integer channel index to extract from x (default 1).
#' @param channels_y Integer channel index to extract from y (default 1).
#' @param spectral_estimator Nonparametric cross-spectral estimator:
#'   `"multitaper"` (default) or `"welch"`.
#' @param n_fft Even FFT length for the nonparametric estimator. `NULL` chooses
#'   a deterministic length from `segment_length`.
#' @param segment_length Segment length for cross-spectral estimation. `NULL`
#'   chooses the largest usable power of two up to 512 samples.
#' @param overlap Fractional segment overlap in `[0, 1)`.
#' @param time_bandwidth DPSS time-bandwidth product for multitaper estimation.
#' @param n_tapers Number of DPSS tapers; `NULL` uses
#'   `floor(2 * time_bandwidth) - 1`.
#' @param factor_tolerance Relative Wilson-factor update tolerance.
#' @param factor_max_iterations Maximum Wilson iterations.
#' @param eigen_floor Relative positive eigenvalue floor used only for
#'   round-off-scale or ill-conditioned cross-spectral bins.
#'
#' @return A list with components:
#'   \describe{
#'     \item{gc_xy}{Numeric Granger causality from x to y (x drives y). A
#'       positive value indicates that x's past helps predict y beyond y's own
#'       past.}
#'     \item{gc_yx}{Numeric Granger causality from y to x (y drives x).}
#'     \item{net_gc}{Numeric net Granger causality (\code{gc_xy - gc_yx}).
#'       Positive values indicate that x drives y more than y drives x.}
#'     \item{order}{Integer AR model order used, or `NA` for nonparametric
#'       estimation.}
#'     \item{method}{Character estimator used.}
#'   }
#'   Nonparametric results additionally contain the frequency grid,
#'   directional spectra, estimator/factorization metadata, innovation
#'   covariance, the factorization reconstruction limit, and regularization
#'   diagnostics.
#'
#' @references
#' Granger, C. W. J. (1969). Investigating causal relations by econometric
#' models and cross-spectral methods. Econometrica, 37(3), 424-438.
#'
#' Geweke, J. (1982). Measurement of linear dependence and feedback between
#' multiple time series. Journal of the American Statistical Association,
#' 77(378), 304-313.
#'
#' Wilson, G. T. (1972). The factorization of matricial spectral densities.
#' SIAM Journal on Applied Mathematics, 23(4), 420-426.
#'
#' Dhamala, M., Rangarajan, G., & Ding, M. (2008). Analyzing information flow
#' in brain networks with nonparametric Granger causality. NeuroImage, 41(2),
#' 354-362.
#'
#' Dhamala, M., Rangarajan, G., & Ding, M. (2008). Estimating Granger
#' causality from Fourier and wavelet transforms of time series data.
#' Physical Review Letters, 100, 018701.
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
                             channels_x = 1L, channels_y = 1L,
                             spectral_estimator = c("multitaper", "welch"),
                             n_fft = NULL, segment_length = NULL,
                             overlap = 0.5, time_bandwidth = 4,
                             n_tapers = NULL, factor_tolerance = 1e-8,
                             factor_max_iterations = 1000L,
                             eigen_floor = 1e-10) {

  method <- match.arg(method)
  if (identical(method, "nonparametric")) {
    pair <- .extract_signal_pair(
      x = x, y = y, sr = sr,
      modality_x = modality_x, modality_y = modality_y,
      channels_x = channels_x, channels_y = channels_y
    )
    return(.gc_nonparametric(
      pair$x, pair$y, pair$sr, freq_range, spectral_estimator,
      n_fft, segment_length, overlap, time_bandwidth, n_tapers,
      factor_tolerance, factor_max_iterations, eigen_floor
    ))
  }

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
    order  = order,
    method = method
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


# ---- MVAR-based DTF / PDC (parity with PhysioEEG eegDTF / eegPDC) -----------

# Assemble a (time x channels) matrix from either a multichannel matrix /
# PhysioExperiment (y = NULL) or a pair of signals (x, y).
.assemble_mvar_input <- function(x, y, sr, modality_x, modality_y,
                                 channels_x, channels_y) {
  if (is.null(y)) {
    if (inherits(x, "PhysioExperiment")) {
      X <- SummarizedExperiment::assay(x, PhysioCore::defaultAssay(x))
      sr <- PhysioCore::samplingRate(x)
      cd <- SummarizedExperiment::colData(x)
      labels <- if ("label" %in% colnames(cd)) as.character(cd$label) else
        paste0("ch", seq_len(ncol(X)))
    } else if (is.matrix(x)) {
      X <- x
      labels <- colnames(x)
      if (is.null(labels)) labels <- paste0("ch", seq_len(ncol(x)))
      if (is.null(sr)) stop("'sr' is required for a matrix input.", call. = FALSE)
    } else {
      stop("Provide 'y', or 'x' as a multichannel matrix / PhysioExperiment.",
           call. = FALSE)
    }
  } else {
    pair <- .extract_signal_pair(
      x, y, sr, modality_x = modality_x, modality_y = modality_y,
      channels_x = if (is.null(channels_x)) 1L else channels_x,
      channels_y = if (is.null(channels_y)) 1L else channels_y)
    X <- cbind(pair$x, pair$y); sr <- pair$sr; labels <- c("x", "y")
  }
  list(X = as.matrix(X), sr = sr, labels = labels)
}

# DTF from H(f): |H_ij| / sqrt(sum_m |H_im|^2) (identical to PhysioEEG eegDTF).
.compute_dtf_cm <- function(H, normalized, ffDTF) {
  n <- dim(H)[1]; nf <- dim(H)[3]; mH <- Mod(H); dtf <- array(0, dim(H))
  if (ffDTF) {
    for (i in seq_len(n)) {
      denom <- sqrt(sum(mH[i, , ]^2))
      dtf[i, , ] <- if (denom > 0) mH[i, , ] / denom else 0
    }
  } else if (normalized) {
    for (f in seq_len(nf)) for (i in seq_len(n)) {
      denom <- sqrt(sum(mH[i, , f]^2))
      dtf[i, , f] <- if (denom > 0) mH[i, , f] / denom else 0
    }
  } else {
    dtf <- mH
  }
  dtf
}

# PDC from Abar(f): |Abar_ij| / sqrt(sum_m |Abar_mj|^2) (identical to eegPDC).
.compute_pdc_cm <- function(Ab, sigma, generalized) {
  n <- dim(Ab)[1]; nf <- dim(Ab)[3]; mA <- Mod(Ab)
  w <- if (generalized) 1 / sigma else rep(1, n); pdc <- array(0, dim(Ab))
  for (f in seq_len(nf)) for (j in seq_len(n)) {
    num <- mA[, j, f] * w; denom <- sqrt(sum(num^2))
    pdc[, j, f] <- if (denom > 0) num / denom else 0
  }
  pdc
}

.band_avg_cm <- function(arr, freqs, band) {
  sel <- if (is.null(band)) rep(TRUE, length(freqs)) else
    freqs >= band[1] & freqs <= band[2]
  if (!any(sel)) sel <- rep(TRUE, length(freqs))
  apply(arr[, , sel, drop = FALSE], c(1, 2), mean)
}

# NA-filled result when the MVAR fit is singular (e.g. a channel against itself
# on the diagonal of a coupling matrix).
.na_directed_result <- function(labels, freqs, n, key, method) {
  bm <- matrix(NA_real_, n, n, dimnames = list(labels, labels))
  arr <- array(NA_real_, c(n, n, length(freqs)),
               dimnames = list(labels, labels, NULL))
  out <- list(matrix = bm, frequencies = freqs, channels = labels,
              order = NA_integer_, method = method)
  out[[key]] <- arr
  if (n == 2L) {
    out[[paste0(key, "_xy")]] <- NA_real_
    out[[paste0(key, "_yx")]] <- NA_real_
  }
  out
}

#' Directed Transfer Function for multichannel / multimodal signals
#'
#' Frequency-resolved directed connectivity via the Directed Transfer Function
#' (Kaminski & Blinowska 1991), computed from the transfer function of a
#' multivariate autoregressive model fitted with [PhysioCore::mvarFit()]. This
#' cross-modal implementation is numerically identical to \code{PhysioEEG}'s
#' \code{eegDTF} on the same input. Because the transfer function inverts the
#' whole VAR system, DTF reflects both direct and indirect pathways.
#'
#' @param x A multichannel matrix (time x channels) or PhysioExperiment (with
#'   \code{y = NULL}), or one signal / PhysioExperiment of a pair.
#' @param y The second signal / PhysioExperiment of a pair, or \code{NULL} for a
#'   single multichannel input.
#' @param sr Sampling rate in Hz (required for numeric-matrix input).
#' @param order MVAR model order (default: 5), or \code{NULL} to select
#'   automatically.
#' @param freqs Frequencies in Hz (default: 128 points from 0 to the Nyquist
#'   frequency).
#' @param normalized Row-normalize the DTF (default: \code{TRUE}).
#' @param ffDTF Use the full-frequency DTF normalization (default: \code{FALSE}).
#' @param band Optional length-2 band in Hz over which to average the returned
#'   matrix (default: all frequencies).
#' @param mvar_method MVAR estimator (default: \code{"ols"}).
#' @param modality_x,modality_y,channels_x,channels_y Passed to the shared
#'   signal-pair extraction when \code{y} is supplied.
#' @return A list with \code{dtf} (n x n x n_freqs array), the band-averaged
#'   directed \code{matrix} (rows = targets, columns = sources),
#'   \code{frequencies}, \code{channels}, and \code{order}. For a two-channel
#'   input, \code{dtf_xy} and \code{dtf_yx} give the directional summaries.
#' @references
#' Kaminski, M. J., & Blinowska, K. J. (1991). A new method of the description
#' of the information flow in the brain structures. Biological Cybernetics,
#' 65(3), 203-210.
#' @seealso [partialDirectedCoherence()], [transferEntropy()], [grangerCausality()]
#' @export
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(1500), ncol = 3)
#' res <- directedTransferFunction(X, sr = 100, order = 4)
#' dim(res$dtf)
directedTransferFunction <- function(x, y = NULL, sr = NULL, order = 5L,
                                     freqs = NULL, normalized = TRUE,
                                     ffDTF = FALSE, band = NULL,
                                     mvar_method = "ols",
                                     modality_x = NULL, modality_y = NULL,
                                     channels_x = NULL, channels_y = NULL) {
  mc <- .assemble_mvar_input(x, y, sr, modality_x, modality_y,
                             channels_x, channels_y)
  Xd <- sweep(mc$X, 2, colMeans(mc$X))
  if (is.null(freqs)) freqs <- seq(0, mc$sr / 2, length.out = 128L)
  fit <- tryCatch(
    PhysioCore::mvarFit(Xd, order = if (is.null(order)) NULL else as.integer(order),
                        method = mvar_method),
    error = function(e) NULL)
  if (is.null(fit)) {                       # singular fit (e.g. collinear input)
    return(.na_directed_result(mc$labels, freqs, ncol(Xd), "dtf",
                               if (ffDTF) "ffDTF" else "dtf"))
  }
  sp <- PhysioCore::mvarTransfer(fit, freqs, mc$sr)
  dtf <- .compute_dtf_cm(sp$H, normalized, ffDTF)
  dimnames(dtf) <- list(mc$labels, mc$labels, NULL)
  bm <- .band_avg_cm(dtf, freqs, band); dimnames(bm) <- list(mc$labels, mc$labels)
  out <- list(dtf = dtf, matrix = bm, frequencies = freqs, channels = mc$labels,
              order = fit$order, method = if (ffDTF) "ffDTF" else "dtf")
  if (ncol(mc$X) == 2L) { out$dtf_xy <- bm[2, 1]; out$dtf_yx <- bm[1, 2] }
  out
}

#' Partial Directed Coherence for multichannel / multimodal signals
#'
#' Frequency-resolved directed connectivity via Partial Directed Coherence
#' (Baccala & Sameshima 2001), computed from the coefficient matrix of an MVAR
#' model ([PhysioCore::mvarFit()]); numerically identical to \code{PhysioEEG}'s
#' \code{eegPDC}. Unlike DTF, PDC reflects only direct influences, so a purely
#' indirect pathway gives PDC near zero.
#'
#' @inheritParams directedTransferFunction
#' @param generalized Use generalized PDC (inverse-noise weighting; default
#'   \code{TRUE}).
#' @return A list as in [directedTransferFunction()] with a \code{pdc} array and
#'   band-averaged \code{matrix}; for a two-channel input, \code{pdc_xy} and
#'   \code{pdc_yx}.
#' @references
#' Baccala, L. A., & Sameshima, K. (2001). Partial directed coherence: a new
#' concept in neural structure determination. Biological Cybernetics, 84(6),
#' 463-474.
#' @seealso [directedTransferFunction()], [transferEntropy()]
#' @export
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(1500), ncol = 3)
#' res <- partialDirectedCoherence(X, sr = 100, order = 4)
#' dim(res$pdc)
partialDirectedCoherence <- function(x, y = NULL, sr = NULL, order = 5L,
                                     freqs = NULL, generalized = TRUE,
                                     band = NULL, mvar_method = "ols",
                                     modality_x = NULL, modality_y = NULL,
                                     channels_x = NULL, channels_y = NULL) {
  mc <- .assemble_mvar_input(x, y, sr, modality_x, modality_y,
                             channels_x, channels_y)
  Xd <- sweep(mc$X, 2, colMeans(mc$X))
  if (is.null(freqs)) freqs <- seq(0, mc$sr / 2, length.out = 128L)
  fit <- tryCatch(
    PhysioCore::mvarFit(Xd, order = if (is.null(order)) NULL else as.integer(order),
                        method = mvar_method),
    error = function(e) NULL)
  if (is.null(fit)) {
    return(.na_directed_result(mc$labels, freqs, ncol(Xd), "pdc",
                               if (generalized) "gpdc" else "pdc"))
  }
  sp <- PhysioCore::mvarTransfer(fit, freqs, mc$sr)
  pdc <- .compute_pdc_cm(sp$A, sqrt(diag(fit$Sigma)), generalized)
  dimnames(pdc) <- list(mc$labels, mc$labels, NULL)
  bm <- .band_avg_cm(pdc, freqs, band); dimnames(bm) <- list(mc$labels, mc$labels)
  out <- list(pdc = pdc, matrix = bm, frequencies = freqs, channels = mc$labels,
              order = fit$order, method = if (generalized) "gpdc" else "pdc")
  if (ncol(mc$X) == 2L) { out$pdc_xy <- bm[2, 1]; out$pdc_yx <- bm[1, 2] }
  out
}


# ---- Transfer entropy internals --------------------------------------------

# Row-wise Chebyshev (max-norm) distance from every row of M to a reference row.
.cheb <- function(M, ref) Reduce(pmax, as.data.frame(abs(sweep(M, 2, ref))))

# Delay embedding: predict tgt[n + delay] from k-history of tgt and l-history of
# src. Returns Yf (future target), Yp (target history), Xp (source history).
.te_embed <- function(src, tgt, k, l, delay) {
  n <- length(tgt); lo <- max(k, l); hi <- n - delay; idx <- lo:hi
  Yp <- vapply(0:(k - 1), function(d) tgt[idx - d], numeric(length(idx)))
  Xp <- vapply(0:(l - 1), function(d) src[idx - d], numeric(length(idx)))
  list(Yf = matrix(tgt[idx + delay], ncol = 1),
       Yp = matrix(Yp, ncol = k), Xp = matrix(Xp, ncol = l))
}

# Histogram (binning) TE = discrete conditional MI I(Yf; Xp | Yp), in bits.
.te_hist <- function(src, tgt, k, l, delay, bins) {
  em <- .te_embed(src, tgt, k, l, delay); n <- nrow(em$Yf)
  code <- function(M) {
    D <- apply(M, 2, function(col) {
      br <- seq(min(col), max(col), length.out = bins + 1)
      as.integer(cut(col, breaks = br, include.lowest = TRUE))
    })
    if (is.null(dim(D))) D <- matrix(D, ncol = 1)
    apply(D, 1, paste, collapse = "_")
  }
  a <- code(em$Yf); b <- code(em$Yp); c <- code(em$Xp)
  pabc <- table(a, b, c) / n; pab <- table(a, b) / n
  pbc <- table(b, c) / n; pb <- table(b) / n
  te <- 0
  for (ia in dimnames(pabc)[[1]]) for (ib in dimnames(pabc)[[2]])
    for (ic in dimnames(pabc)[[3]]) {
      p <- pabc[ia, ib, ic]; if (p <= 0) next
      den <- pab[ia, ib] * pbc[ib, ic]; if (den <= 0) next
      te <- te + p * log2(p * pb[ib] / den)
    }
  as.numeric(te)
}

# KSG (Kraskov) k-NN TE = continuous conditional MI I(Yf; Xp | Yp), in nats.
.te_ksg <- function(src, tgt, k, l, delay, knn) {
  em <- .te_embed(src, tgt, k, l, delay)
  Yf <- em$Yf; Yp <- em$Yp; Xp <- em$Xp; N <- nrow(Yf)
  W <- cbind(Yf, Yp, Xp); YfYp <- cbind(Yf, Yp); XpYp <- cbind(Xp, Yp)
  acc <- numeric(N)
  for (i in seq_len(N)) {
    eps <- sort(.cheb(W, W[i, ]))[knn + 1]
    n_yp <- sum(.cheb(Yp, Yp[i, ]) < eps) - 1
    n_yfyp <- sum(.cheb(YfYp, YfYp[i, ]) < eps) - 1
    n_xpyp <- sum(.cheb(XpYp, XpYp[i, ]) < eps) - 1
    acc[i] <- digamma(n_yp + 1) - digamma(n_yfyp + 1) - digamma(n_xpyp + 1)
  }
  digamma(knn) + mean(acc)
}

#' Transfer entropy (model-free directed information)
#'
#' Estimates the transfer entropy from \code{x} to \code{y} and back
#' (Schreiber 2000), a model-free measure of directed information flow, with a
#' histogram (binning) estimator or the Kraskov-Stoegbauer-Grassberger (KSG)
#' k-nearest-neighbour estimator. With \code{effective = TRUE} the estimate is
#' bias-corrected by subtracting the mean transfer entropy of IAAFT surrogates
#' of the source, so independent signals give an effective transfer entropy near
#' zero.
#'
#' @param x,y Two numeric signals or PhysioExperiment objects (or a
#'   MultiPhysioExperiment via \code{x}).
#' @param sr Sampling rate in Hz (required for numeric input).
#' @param k Target (destination) embedding length (default: 1).
#' @param l Source embedding length (default: 1).
#' @param delay Source-to-target lag in samples (default: 1).
#' @param estimator \code{"ksg"} (k-NN, in nats) or \code{"histogram"}
#'   (binning, in bits).
#' @param knn Number of neighbours for the KSG estimator (default: 4).
#' @param bins Number of bins per dimension for the histogram estimator
#'   (default: 8).
#' @param effective Subtract the surrogate (IAAFT) null mean (default:
#'   \code{FALSE}).
#' @param n_surrogate Number of surrogates for \code{effective} (default: 19).
#' @param seed Optional RNG seed for reproducible surrogates.
#' @param modality_x,modality_y,channels_x,channels_y Passed to the shared
#'   signal-pair extraction.
#' @return A list with \code{te_xy}, \code{te_yx}, and \code{net}
#'   (\code{te_xy - te_yx}); when \code{effective = TRUE} also \code{eff_xy},
#'   \code{eff_yx}, and the surrogate null standard deviations.
#' @references
#' Schreiber, T. (2000). Measuring information transfer. Physical Review
#' Letters, 85(2), 461-464.
#'
#' Kraskov, A., Stoegbauer, H., & Grassberger, P. (2004). Estimating mutual
#' information. Physical Review E, 69(6), 066138.
#' @seealso [directedTransferFunction()], [grangerCausality()],
#'   [surrogateTest()]
#' @export
#' @examples
#' set.seed(1)
#' x <- rnorm(600); y <- c(rep(0, 5), 0.7 * x[1:595]) + rnorm(600, sd = 0.5)
#' transferEntropy(x, y, sr = 100, delay = 5, estimator = "histogram")$net
transferEntropy <- function(x, y = NULL, sr = NULL, k = 1L, l = 1L, delay = 1L,
                            estimator = c("ksg", "histogram"), knn = 4L,
                            bins = 8L, effective = FALSE, n_surrogate = 19L,
                            seed = NULL, modality_x = NULL, modality_y = NULL,
                            channels_x = 1L, channels_y = 1L) {
  estimator <- match.arg(estimator)
  pair <- .extract_signal_pair(x, y, sr, modality_x = modality_x,
                               modality_y = modality_y,
                               channels_x = channels_x, channels_y = channels_y)
  sx <- pair$x; sy <- pair$y
  est <- if (estimator == "ksg") {
    function(a, b) .te_ksg(a, b, k, l, delay, knn)
  } else {
    function(a, b) .te_hist(a, b, k, l, delay, bins)
  }
  te_xy <- est(sx, sy); te_yx <- est(sy, sx)
  out <- list(te_xy = te_xy, te_yx = te_yx, net = te_xy - te_yx,
              estimator = estimator, k = k, l = l, delay = delay)
  if (effective) {
    if (!is.null(seed)) set.seed(seed)
    sxy <- vapply(seq_len(n_surrogate),
                  function(i) est(.iaaft_surrogate(sx), sy), numeric(1))
    syx <- vapply(seq_len(n_surrogate),
                  function(i) est(.iaaft_surrogate(sy), sx), numeric(1))
    out$eff_xy <- te_xy - mean(sxy); out$eff_yx <- te_yx - mean(syx)
    out$null_sd_xy <- stats::sd(sxy); out$null_sd_yx <- stats::sd(syx)
  }
  out
}
