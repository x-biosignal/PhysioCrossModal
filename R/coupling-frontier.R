# Frontier connectivity measures ------------------------------------------
# Volume-conduction / leakage-robust and directed measures that complement the
# coherence / PLV / directed-coupling suite: the Phase Slope Index (directed,
# VC-robust), orthogonalised amplitude-envelope correlation (leakage-corrected),
# corrected imaginary PLV, pairwise phase consistency, and LEiDA dynamic
# functional-connectivity states.

#' @keywords internal
#' @noRd
.frontier_pair <- function(x, y, sr, modality_x, modality_y,
                           channels_x, channels_y) {
  pair <- .extract_signal_pair(x, y, sr = sr,
                               modality_x = modality_x, modality_y = modality_y,
                               channels_x = channels_x, channels_y = channels_y)
  n <- min(length(pair$x), length(pair$y))
  list(x = pair$x[seq_len(n)], y = pair$y[seq_len(n)], sr = pair$sr)
}

#' @keywords internal
#' @noRd
.analytic_band <- function(sig, band, sr, order = 4L) {
  .hilbert_analytic(.bandpass_filter(sig, band[1], band[2], sr, order = order))
}

#' @keywords internal
#' @noRd
.check_band2 <- function(freq_band) {
  if (missing(freq_band) || !is.numeric(freq_band) || length(freq_band) != 2 ||
      freq_band[1] >= freq_band[2])
    stop("'freq_band' must be a numeric c(low, high) with low < high.",
         call. = FALSE)
}

#' Phase Slope Index (directed, volume-conduction robust)
#'
#' The Phase Slope Index (Nolte et al. 2008) estimates the direction of
#' information flow from the slope of the coherency phase across a frequency
#' band, and is insensitive to instantaneous (zero-lag) mixing. Positive values
#' indicate that `x` leads (drives) `y`.
#'
#' @param x,y Signals (or a single object from which a pair is extracted).
#' @param sr Sampling rate in Hz.
#' @param freq_band Numeric `c(low, high)` band in Hz.
#' @param nperseg,noverlap Welch segment length / overlap for the cross-spectrum.
#' @param normalize If `TRUE` (default), scale into `[-1, 1]` by the summed
#'   coherency magnitude; if `FALSE`, return the raw imaginary sum.
#' @param modality_x,modality_y,channels_x,channels_y Passed to the internal
#'   signal-pair extractor.
#' @return A single Phase Slope Index value.
#' @references Nolte G, et al. (2008). Robustly estimating the flow direction of
#'   information in complex physical systems. \emph{Phys Rev Lett} 100:234101.
#'   \doi{10.1103/PhysRevLett.100.234101}
#' @examples
#' set.seed(1); sr <- 200; n <- 4000
#' x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2)); x[is.na(x)] <- 0
#' y <- c(rep(0, 8), x[seq_len(n - 8)])          # y lags x
#' phaseSlopeIndex(x, y, sr = sr, freq_band = c(1, 40))
#' @export
phaseSlopeIndex <- function(x, y = NULL, sr = NULL, freq_band,
                            nperseg = 256L, noverlap = NULL, normalize = TRUE,
                            modality_x = NULL, modality_y = NULL,
                            channels_x = 1L, channels_y = 1L) {
  .check_band2(freq_band)
  p <- .frontier_pair(x, y, sr, modality_x, modality_y, channels_x, channels_y)
  Sxy <- .welch_csd(p$x, p$y, nperseg = nperseg, noverlap = noverlap, sr = p$sr)
  Sxx <- Re(.welch_csd(p$x, p$x, nperseg = nperseg, noverlap = noverlap,
                       sr = p$sr)$csd)
  Syy <- Re(.welch_csd(p$y, p$y, nperseg = nperseg, noverlap = noverlap,
                       sr = p$sr)$csd)
  f <- Sxy$frequencies
  C <- Sxy$csd / sqrt(Sxx * Syy)                 # complex coherency
  idx <- which(f >= freq_band[1] & f <= freq_band[2])
  if (length(idx) < 2)
    stop("'freq_band' spans too few frequency bins; widen it or raise nperseg.",
         call. = FALSE)
  a <- idx[-length(idx)]; b <- idx[-1]
  psi <- Im(sum(Conj(C[a]) * C[b]))
  if (normalize) {
    den <- sum(Mod(C[a]) * Mod(C[b]))
    psi <- if (den > 0) psi / den else 0
  }
  psi
}

#' Orthogonalised amplitude-envelope correlation (leakage-corrected)
#'
#' Amplitude-envelope correlation after pairwise orthogonalisation (Hipp et al.
#' 2012; Brookes et al. 2012): the instantaneous zero-lag (leakage) component is
#' projected out before correlating band-limited amplitude envelopes, so shared
#' instantaneous mixing does not create spurious coupling. Symmetrised over the
#' two orthogonalisation directions.
#'
#' @inheritParams phaseSlopeIndex
#' @param order Butterworth filter order for the band-pass (default 4).
#' @return A single leakage-corrected envelope correlation in `[-1, 1]`.
#' @references Hipp JF, et al. (2012). Large-scale cortical correlation structure
#'   of spontaneous oscillatory activity. \emph{Nat Neurosci} 15:884-890.
#' @examples
#' set.seed(1); sr <- 200; n <- 4000
#' x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2)); x[is.na(x)] <- 0
#' y <- c(rep(0, 10), x[seq_len(n - 10)]) + rnorm(n, sd = 0.5)
#' orthogonalizedAEC(x, y, sr = sr, freq_band = c(8, 12))
#' @export
orthogonalizedAEC <- function(x, y = NULL, sr = NULL, freq_band, order = 4L,
                              modality_x = NULL, modality_y = NULL,
                              channels_x = 1L, channels_y = 1L) {
  .check_band2(freq_band)
  p <- .frontier_pair(x, y, sr, modality_x, modality_y, channels_x, channels_y)
  ax <- .analytic_band(p$x, freq_band, p$sr, order)
  ay <- .analytic_band(p$y, freq_band, p$sr, order)
  orth_env <- function(a, ref) {
    proj <- Re(a * Conj(ref) / (Mod(ref)^2 + .Machine$double.eps))
    Mod(a - proj * ref)
  }
  safe_cor <- function(a, b) {
    if (stats::sd(a) == 0 || stats::sd(b) == 0) return(0)
    stats::cor(a, b)
  }
  aec1 <- safe_cor(Mod(ax), orth_env(ay, ax))
  aec2 <- safe_cor(Mod(ay), orth_env(ax, ay))
  (aec1 + aec2) / 2
}

#' Corrected imaginary phase-locking value (ciPLV)
#'
#' ciPLV (Bruna et al. 2018) corrects the imaginary part of the PLV for the size
#' of its real part, giving a phase-synchrony estimate that is robust to zero-lag
#' (volume-conduction) coupling.
#'
#' @inheritParams orthogonalizedAEC
#' @return A single ciPLV magnitude in `[0, 1]`.
#' @references Bruna R, et al. (2018). Phase locking value revisited. \emph{J
#'   Neural Eng} 15:056011.
#' @export
ciPLV <- function(x, y = NULL, sr = NULL, freq_band, order = 4L,
                  modality_x = NULL, modality_y = NULL,
                  channels_x = 1L, channels_y = 1L) {
  .check_band2(freq_band)
  p <- .frontier_pair(x, y, sr, modality_x, modality_y, channels_x, channels_y)
  ax <- .analytic_band(p$x, freq_band, p$sr, order)
  ay <- .analytic_band(p$y, freq_band, p$sr, order)
  z <- mean(exp(1i * (Arg(ax) - Arg(ay))))
  denom <- sqrt(1 - Re(z)^2)
  if (denom <= 0) return(0)
  abs(Im(z) / denom)
}

#' Pairwise phase consistency (PPC)
#'
#' The bias-free phase-synchrony estimator of Vinck et al. (2010):
#' \eqn{PPC = (|\sum e^{i\Delta\phi}|^2 - N) / (N(N-1))}. Unlike the PLV it is not
#' inflated by the sample size.
#'
#' @inheritParams orthogonalizedAEC
#' @return A single PPC value (approaches 1 for perfect locking, 0 for none).
#' @references Vinck M, et al. (2010). The pairwise phase consistency.
#'   \emph{NeuroImage} 51:112-122.
#' @export
pairwisePhaseConsistency <- function(x, y = NULL, sr = NULL, freq_band,
                                     order = 4L,
                                     modality_x = NULL, modality_y = NULL,
                                     channels_x = 1L, channels_y = 1L) {
  .check_band2(freq_band)
  p <- .frontier_pair(x, y, sr, modality_x, modality_y, channels_x, channels_y)
  ax <- .analytic_band(p$x, freq_band, p$sr, order)
  ay <- .analytic_band(p$y, freq_band, p$sr, order)
  z <- exp(1i * (Arg(ax) - Arg(ay)))
  N <- length(z)
  (Mod(sum(z))^2 - N) / (N * (N - 1))
}

#' LEiDA dynamic functional-connectivity states
#'
#' Leading Eigenvector Dynamics Analysis (Cabral et al. 2017): at each time
#' point the instantaneous phase-coherence matrix is reduced to its leading
#' eigenvector; these are clustered over time into a small set of recurrent
#' connectivity states. Complements state-resolved network analysis by finding
#' the states data-drivenly rather than from external labels.
#'
#' @param x A numeric matrix, time (rows) by channels (columns).
#' @param sr Sampling rate in Hz.
#' @param freq_band Numeric `c(low, high)` band in Hz.
#' @param n_states Number of clusters/states (default 4).
#' @param order Butterworth filter order (default 4).
#' @param nstart k-means restarts (default 10).
#' @param seed Optional integer seed for reproducible clustering.
#' @return An object of class `"leida"`: a list with `states` (per-time cluster
#'   labels), `centroids` (states x channels leading-eigenvector centroids),
#'   `occupancy` (fraction of time in each state), `dwell` (mean run length in
#'   samples per state), and `n_states`.
#' @references Cabral J, et al. (2017). Cognitive performance in healthy older
#'   adults relates to spontaneous switching between states of functional
#'   connectivity during rest. \emph{Sci Rep} 7:5135.
#' @examples
#' set.seed(1); sr <- 100; n <- 2000; t <- seq_len(n) / sr
#' # two regimes: channels 1-2 coherent first half, 3-4 coherent second half
#' base <- sin(2 * pi * 10 * t)
#' X <- cbind(base + rnorm(n, sd = .3), base + rnorm(n, sd = .3),
#'            sin(2 * pi * 10 * t + 1) + rnorm(n, sd = .3),
#'            sin(2 * pi * 10 * t + 1) + rnorm(n, sd = .3))
#' res <- leidaStates(X, sr = sr, freq_band = c(8, 12), n_states = 2, seed = 1)
#' res$occupancy
#' @export
leidaStates <- function(x, sr, freq_band, n_states = 4L, order = 4L,
                        nstart = 10L, seed = NULL) {
  .check_band2(freq_band)
  X <- as.matrix(x)
  if (ncol(X) < 2L) stop("`x` must have at least two channels (columns).",
                         call. = FALSE)
  nt <- nrow(X); nch <- ncol(X)
  ph <- apply(X, 2L, function(s)
    Arg(.analytic_band(s, freq_band, sr, order)))
  V <- matrix(0, nt, nch)
  for (t in seq_len(nt)) {
    dPC <- cos(outer(ph[t, ], ph[t, ], "-"))
    v1 <- eigen(dPC, symmetric = TRUE)$vectors[, 1]
    if (sum(v1 > 0) < sum(v1 < 0)) v1 <- -v1        # sign convention
    V[t, ] <- v1
  }
  if (!is.null(seed)) set.seed(seed)
  km <- stats::kmeans(V, centers = n_states, nstart = nstart)
  states <- km$cluster
  r <- rle(states)
  dwell <- vapply(seq_len(n_states), function(s) {
    runs <- r$lengths[r$values == s]
    if (length(runs)) mean(runs) else 0
  }, numeric(1))
  structure(
    list(states = states, centroids = km$centers,
         occupancy = tabulate(states, n_states) / nt,
         dwell = dwell, n_states = n_states),
    class = "leida")
}

#' @export
print.leida <- function(x, ...) {
  cat("<LEiDA dynamic connectivity states>\n")
  cat(sprintf("  states:    %d over %d samples\n",
              x$n_states, length(x$states)))
  cat(sprintf("  occupancy: %s\n",
              paste(sprintf("%.2f", x$occupancy), collapse = " ")))
  cat(sprintf("  mean dwell: %s samples\n",
              paste(sprintf("%.0f", x$dwell), collapse = " ")))
  invisible(x)
}
