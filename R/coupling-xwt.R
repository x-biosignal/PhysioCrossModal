# Standalone cross-wavelet transform with red-noise significance ----------
# waveletCoherence() exposes coherence; this adds the cross-wavelet transform
# (XWT) as a first-class result -- common oscillatory power and phase lead/lag
# between two signals in time-frequency -- with an AR(1) red-noise significance
# test (Torrence & Compo 1998; Grinsted et al. 2004). Significance is assessed
# by Monte-Carlo against AR(1) surrogates pushed through the same wavelet
# transform, so it is independent of the internal wavelet normalisation.

#' @keywords internal
#' @noRd
.lag1_acf <- function(v) {
  a <- tryCatch(stats::acf(v, lag.max = 1, plot = FALSE)$acf[2],
                error = function(e) 0)
  if (!is.finite(a)) 0 else max(min(a, 0.999), 0)
}

# AR(1) red-noise surrogate with a given lag-1 autocorrelation and sd.
#' @keywords internal
#' @noRd
.ar1_surrogate <- function(n, alpha, sdv) {
  e <- stats::rnorm(n)
  s <- numeric(n); s[1] <- e[1]
  for (i in seq_len(n)[-1]) s[i] <- alpha * s[i - 1] + e[i]
  s <- s - mean(s)
  sds <- stats::sd(s)
  if (sds > 0) s <- s / sds * sdv
  s
}

#' Cross-wavelet transform (XWT) with red-noise significance
#'
#' Computes the cross-wavelet transform \eqn{W_{xy} = W_x W_y^*} of two signals
#' via the Morlet wavelet: its modulus is the common time-frequency power and its
#' argument is the phase lead/lag between the signals. Common power is tested
#' against an AR(1) red-noise background (Torrence & Compo 1998; Grinsted et al.
#' 2004) by Monte-Carlo: AR(1) surrogates matched to each signal's lag-1
#' autocorrelation and variance are transformed identically, and the per-
#' frequency `siglvl` quantile of the surrogate cross-power sets the threshold.
#'
#' @param x,y Signals (or a single object from which a pair is extracted).
#' @param sr Sampling rate in Hz.
#' @param frequencies Frequencies (Hz) to analyse (default `seq(1, 40, 1)`).
#' @param n_cycles Morlet wavelet cycles (default 7).
#' @param significance If `TRUE` (default), compute red-noise significance.
#' @param n_surrogates Number of AR(1) surrogate pairs (default 100).
#' @param siglvl Significance quantile (default 0.95).
#' @param modality_x,modality_y,channels_x,channels_y Passed to the pair extractor.
#'
#' @return An object of class `"cross_wavelet"`: a list with `power`
#'   (`time x frequency` cross-wavelet modulus), `phase` (relative phase, radians;
#'   positive = `x` leads `y`), `cross` (complex XWT), `frequencies`, `times`,
#'   `coi` (cone of influence), and — when `significance` — `sig_level` (per
#'   frequency), `sig_ratio` (power / level), and `significant` (logical mask).
#'
#' @references
#' Torrence C, Compo GP (1998). A practical guide to wavelet analysis.
#' \emph{Bull Amer Meteor Soc} 79:61-78.
#' Grinsted A, Moore JC, Jevrejeva S (2004). Application of the cross wavelet
#' transform and wavelet coherence to geophysical time series.
#' \emph{Nonlin Processes Geophys} 11:561-566.
#'
#' @examples
#' set.seed(1); sr <- 100; t <- seq(0, 20, by = 1 / sr)
#' x <- sin(2 * pi * 10 * t) + rnorm(length(t), sd = 0.5)
#' y <- sin(2 * pi * 10 * t - pi / 2) + rnorm(length(t), sd = 0.5)  # y lags x
#' xw <- crossWaveletTransform(x, y, sr = sr, frequencies = seq(4, 16, by = 1),
#'                             n_surrogates = 50)
#' mean(xw$significant)
#' @export
crossWaveletTransform <- function(x, y = NULL, sr = NULL,
                                  frequencies = seq(1, 40, by = 1),
                                  n_cycles = 7, significance = TRUE,
                                  n_surrogates = 100L, siglvl = 0.95,
                                  modality_x = NULL, modality_y = NULL,
                                  channels_x = 1L, channels_y = 1L) {
  pair <- .extract_signal_pair(x, y, sr = sr,
                               modality_x = modality_x, modality_y = modality_y,
                               channels_x = channels_x, channels_y = channels_y)
  n <- min(length(pair$x), length(pair$y))
  sx <- pair$x[seq_len(n)]; sy <- pair$y[seq_len(n)]
  srv <- pair$sr

  W_x <- .morlet_wavelet_transform(sx, frequencies, n_cycles, srv)
  W_y <- .morlet_wavelet_transform(sy, frequencies, n_cycles, srv)
  W_xy <- W_x * Conj(W_y)

  out <- list(power = Mod(W_xy), phase = Arg(W_xy), cross = W_xy,
              frequencies = frequencies, times = seq_len(n) / srv,
              coi = .compute_coi(n, srv, n_cycles))

  if (significance) {
    ax <- .lag1_acf(sx); ay <- .lag1_acf(sy)
    sdx <- stats::sd(sx); sdy <- stats::sd(sy)
    nf <- length(frequencies)
    sig_q <- matrix(0, n_surrogates, nf)
    for (m in seq_len(n_surrogates)) {
      xs <- .ar1_surrogate(n, ax, sdx)
      ys <- .ar1_surrogate(n, ay, sdy)
      Ps <- Mod(.morlet_wavelet_transform(xs, frequencies, n_cycles, srv) *
                Conj(.morlet_wavelet_transform(ys, frequencies, n_cycles, srv)))
      sig_q[m, ] <- apply(Ps, 2L, stats::quantile, probs = siglvl, names = FALSE)
    }
    sig_level <- colMeans(sig_q)                       # per frequency
    sigmat <- matrix(sig_level, nrow = n, ncol = nf, byrow = TRUE)
    out$sig_level <- sig_level
    out$sig_ratio <- out$power / sigmat
    out$significant <- out$sig_ratio > 1
  }
  structure(out, class = "cross_wavelet")
}

#' @export
print.cross_wavelet <- function(x, ...) {
  cat("<Cross-wavelet transform>\n")
  cat(sprintf("  time-frequency: %d samples x %d frequencies (%.1f-%.1f Hz)\n",
              nrow(x$power), length(x$frequencies),
              min(x$frequencies), max(x$frequencies)))
  pk <- which(x$power == max(x$power), arr.ind = TRUE)[1, ]
  cat(sprintf("  peak power at %.1f Hz, t = %.2f s\n",
              x$frequencies[pk[2]], x$times[pk[1]]))
  if (!is.null(x$significant))
    cat(sprintf("  significant (vs red noise): %.1f%% of the plane\n",
                100 * mean(x$significant)))
  invisible(x)
}
