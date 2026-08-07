# Wavelet-based coupling analysis
#
# Implements time-frequency coherence and PLV using Morlet wavelet
# transforms, providing localised coupling measures in both time and
# frequency.

# ---- Morlet wavelet transform ------------------------------------------------

#' Complex Morlet wavelet transform via FFT convolution
#'
#' Computes the continuous wavelet transform of a signal using complex
#' Morlet wavelets, one per requested frequency. Convolution is performed
#' in the frequency domain (FFT) for efficiency.
#'
#' @param x Numeric vector — input signal.
#' @param frequencies Numeric vector of centre frequencies (Hz).
#' @param n_cycles Numeric number of wavelet cycles (higher = better
#'   frequency resolution, lower = better temporal resolution). Default 7.
#' @param sr Numeric sampling rate in Hz.
#' @return Complex matrix of dimension \code{\code{[time x frequency]}}.
#' @noRd
.morlet_wavelet_transform <- function(x, frequencies, n_cycles, sr) {
  n <- length(x)
  n_freq <- length(frequencies)

  # FFT of the signal (zero-pad for linear convolution)
  n_fft <- stats::nextn(n + max(n, 1024L), factors = 2L)
  fft_x <- stats::fft(c(x, rep(0, n_fft - n)))

  # Frequency axis for FFT bins
  fft_freqs <- seq(0, sr * (1 - 1 / n_fft), length.out = n_fft)

  result <- matrix(complex(0), nrow = n, ncol = n_freq)

  for (fi in seq_len(n_freq)) {
    f0 <- frequencies[fi]
    if (f0 <= 0) next

    # Morlet wavelet in frequency domain:
    #   W(f) = exp(-0.5 * ((f - f0) / sigma_f)^2)
    # where sigma_f = f0 / n_cycles
    sigma_f <- f0 / n_cycles
    gauss <- exp(-0.5 * ((fft_freqs - f0) / sigma_f)^2)
    # Zero out negative frequency contributions for analytic signal
    # (fft_freqs > sr/2 correspond to negative frequencies)
    gauss[fft_freqs > sr / 2] <- 0

    # Convolve via multiplication in frequency domain
    conv <- stats::fft(fft_x * gauss, inverse = TRUE) / n_fft
    result[, fi] <- conv[seq_len(n)]
  }

  result
}

# ---- Cone of Influence -------------------------------------------------------

#' Compute Cone of Influence for Morlet wavelet
#'
#' The COI defines the region where edge effects are significant. For Morlet
#' wavelets, the e-folding time is \code{sqrt(2) * n_cycles / (2 * pi * f)},
#' so the COI at each time point is the minimum frequency resolvable without
#' edge contamination.
#'
#' @param n Integer signal length.
#' @param sr Numeric sampling rate.
#' @param n_cycles Numeric wavelet cycles.
#' @return Numeric vector of length \code{n}: the COI frequency at each time
#'   point. Frequencies below this curve are affected by edge artifacts.
#' @noRd
.compute_coi <- function(n, sr, n_cycles) {
  t <- seq(0, (n - 1) / sr, length.out = n)
  # Distance to nearest edge (in seconds)
  edge_dist <- pmin(t, t[n] - t)
  # COI frequency: f = sqrt(2) * n_cycles / (2 * pi * edge_dist)
  # At the edges, edge_dist -> 0, so COI freq -> Inf (everything unreliable)
  coi_freq <- sqrt(2) * n_cycles / (2 * pi * pmax(edge_dist, 1 / sr))
  coi_freq
}

# ---- Exported functions ------------------------------------------------------

#' Time-frequency wavelet coherence
#'
#' Computes wavelet coherence between two signals using complex Morlet
#' wavelets. The cross-wavelet spectrum and auto-spectra are smoothed
#' with a Gaussian temporal window (width proportional to
#' \code{smoothing_cycles / frequency}), and coherence is computed as:
#'
#' \deqn{C(t,f) = \frac{|\langle W_{xy}(t,f) \rangle|^2}{\langle |W_x(t,f)|^2 \rangle \cdot \langle |W_y(t,f)|^2 \rangle}}
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when \code{x} is an
#'   MPE.
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param frequencies Numeric vector of centre frequencies in Hz
#'   (default \code{seq(1, 40, by = 1)}).
#' @param n_cycles Numeric number of wavelet cycles (default 7).
#' @param smoothing_cycles Numeric number of cycles for the temporal
#'   smoothing Gaussian (default 3).
#' @param modality_x,modality_y Character modality names for MPE input.
#' @param channels_x,channels_y Integer channel indices (default 1).
#' @param ... Currently unused.
#'
#' @return A list with components:
#'   \describe{
#'     \item{coherence}{Numeric matrix \code{[time x frequency]} of coherence
#'       values in \eqn{[0, 1]}.}
#'     \item{phase}{Numeric matrix \code{[time x frequency]} of phase differences
#'       (radians).}
#'     \item{frequencies}{Numeric vector of centre frequencies.}
#'     \item{times}{Numeric vector of time points (seconds from start).}
#'     \item{coi}{Numeric vector of Cone of Influence frequencies. At each
#'       time point, frequencies below this value are affected by edge
#'       artifacts.}
#'   }
#'
#' @references
#' Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet
#' analysis. \emph{Bulletin of the American Meteorological Society},
#' 79(1), 61--78.
#'
#' Grinsted, A., Moore, J. C., & Jevrejeva, S. (2004). Application of the
#' cross wavelet transform and wavelet coherence to geophysical time series.
#' \emph{Nonlinear Processes in Geophysics}, 11(5/6), 561--566.
#'
#' @seealso [waveletPLV()], [coherence()], [plotWaveletCoherence()],
#'   [couplingAnalysis()]
#' @export
#' @examples
#' sr <- 200
#' t <- seq(0, 2, length.out = sr * 2)
#' x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
#' y <- 0.8 * sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
#' result <- waveletCoherence(x, y, sr = sr, frequencies = seq(5, 20))
waveletCoherence <- function(x, y = NULL, sr = NULL,
                             frequencies = seq(1, 40, by = 1),
                             n_cycles = 7, smoothing_cycles = 3,
                             modality_x = NULL, modality_y = NULL,
                             channels_x = 1L, channels_y = 1L,
                             ...) {

  # Extract signal pair

  pair <- .extract_signal_pair(x, y, sr = sr,
                               modality_x = modality_x,
                               modality_y = modality_y,
                               channels_x = channels_x,
                               channels_y = channels_y)
  sig_x <- pair$x
  sig_y <- pair$y
  sr_val <- pair$sr

  # Ensure equal length
  n <- min(length(sig_x), length(sig_y))
  sig_x <- sig_x[seq_len(n)]
  sig_y <- sig_y[seq_len(n)]

  # Compute wavelet transforms
  W_x <- .morlet_wavelet_transform(sig_x, frequencies, n_cycles, sr_val)
  W_y <- .morlet_wavelet_transform(sig_y, frequencies, n_cycles, sr_val)

  # Cross- and auto-spectra
  W_xy <- W_x * Conj(W_y)
  P_xx <- Mod(W_x)^2
  P_yy <- Mod(W_y)^2

  n_freq <- length(frequencies)
  coh <- matrix(0, nrow = n, ncol = n_freq)
  phase_diff <- matrix(0, nrow = n, ncol = n_freq)

  for (fi in seq_len(n_freq)) {
    f0 <- frequencies[fi]
    # Gaussian smoothing width in samples: smoothing_cycles / f0 * sr
    sigma_t <- smoothing_cycles / f0 * sr_val
    # Truncate kernel at 3 sigma
    half_width <- min(ceiling(3 * sigma_t), floor(n / 2))
    kernel_idx <- seq(-half_width, half_width)
    kernel <- exp(-0.5 * (kernel_idx / sigma_t)^2)
    kernel <- kernel / sum(kernel)

    # Smooth cross-spectrum and auto-spectra
    smooth_xy <- .convolve_mirror(W_xy[, fi], kernel)
    smooth_xx <- .convolve_mirror(P_xx[, fi], kernel)
    smooth_yy <- .convolve_mirror(P_yy[, fi], kernel)

    denom <- smooth_xx * smooth_yy
    valid <- denom > .Machine$double.eps
    coh[valid, fi] <- Mod(smooth_xy[valid])^2 / denom[valid]
    phase_diff[, fi] <- Arg(smooth_xy)
  }

  # Clamp to \eqn{[0, 1]}
  coh <- pmin(pmax(coh, 0), 1)

  times <- seq(0, (n - 1) / sr_val, length.out = n)

  list(
    coherence = coh,
    phase = phase_diff,
    frequencies = frequencies,
    times = times,
    coi = .compute_coi(n, sr_val, n_cycles)
  )
}


#' Time-frequency wavelet PLV
#'
#' Computes Phase Locking Value in the time-frequency domain using
#' complex Morlet wavelets. The phase difference between the two signals
#' is computed at each time-frequency point, and PLV is the magnitude of
#' the smoothed unit-phase vector:
#'
#' \deqn{\text{PLV}(t,f) = \left|\langle e^{i\Delta\phi(t,f)} \rangle\right|}
#'
#' @inheritParams waveletCoherence
#'
#' @return A list with components:
#'   \describe{
#'     \item{plv}{Numeric matrix \code{[time x frequency]} of PLV values
#'       in \eqn{[0, 1]}.}
#'     \item{frequencies}{Numeric vector of centre frequencies.}
#'     \item{times}{Numeric vector of time points (seconds from start).}
#'     \item{coi}{Numeric vector of Cone of Influence frequencies. At each
#'       time point, frequencies below this value are affected by edge
#'       artifacts.}
#'   }
#'
#' @references
#' Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet
#' analysis. \emph{Bulletin of the American Meteorological Society},
#' 79(1), 61--78.
#'
#' Lachaux, J.-P., Rodriguez, E., Martinerie, J., & Varela, F. J. (1999).
#' Measuring phase synchrony in brain signals. \emph{Human Brain Mapping},
#' 8(4), 194--208.
#'
#' @seealso [waveletCoherence()], [phaseLockingValue()],
#'   [plotWaveletCoherence()], [couplingAnalysis()]
#' @export
#' @examples
#' sr <- 200
#' t <- seq(0, 2, length.out = sr * 2)
#' x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
#' y <- sin(2 * pi * 10 * t + pi/4) + 0.3 * rnorm(length(t))
#' result <- waveletPLV(x, y, sr = sr, frequencies = seq(5, 20))
waveletPLV <- function(x, y = NULL, sr = NULL,
                       frequencies = seq(1, 40, by = 1),
                       n_cycles = 7, smoothing_cycles = 3,
                       modality_x = NULL, modality_y = NULL,
                       channels_x = 1L, channels_y = 1L,
                       ...) {

  # Extract signal pair
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

  # Compute wavelet transforms
  W_x <- .morlet_wavelet_transform(sig_x, frequencies, n_cycles, sr_val)
  W_y <- .morlet_wavelet_transform(sig_y, frequencies, n_cycles, sr_val)

  # Phase difference unit vectors
  phase_diff_vec <- W_x * Conj(W_y)
  # Normalize to unit magnitude
  magnitudes <- Mod(phase_diff_vec)
  magnitudes[magnitudes < .Machine$double.eps] <- 1
  unit_vectors <- phase_diff_vec / magnitudes

  n_freq <- length(frequencies)
  plv_mat <- matrix(0, nrow = n, ncol = n_freq)

  for (fi in seq_len(n_freq)) {
    f0 <- frequencies[fi]
    sigma_t <- smoothing_cycles / f0 * sr_val
    half_width <- min(ceiling(3 * sigma_t), floor(n / 2))
    kernel_idx <- seq(-half_width, half_width)
    kernel <- exp(-0.5 * (kernel_idx / sigma_t)^2)
    kernel <- kernel / sum(kernel)

    smoothed <- .convolve_mirror(unit_vectors[, fi], kernel)
    plv_mat[, fi] <- Mod(smoothed)
  }

  # Clamp to \eqn{[0, 1]}
  plv_mat <- pmin(pmax(plv_mat, 0), 1)

  times <- seq(0, (n - 1) / sr_val, length.out = n)

  list(
    plv = plv_mat,
    frequencies = frequencies,
    times = times,
    coi = .compute_coi(n, sr_val, n_cycles)
  )
}


# ---- Internal smoothing helper -----------------------------------------------

#' Convolve with mirror-padded boundaries
#'
#' Convolves a (possibly complex) vector with a symmetric kernel, using
#' mirror-reflection at the edges to reduce boundary artifacts.
#'
#' @param x Numeric or complex vector.
#' @param kernel Numeric vector (symmetric smoothing kernel).
#' @return Vector of the same length and type as \code{x}.
#' @noRd
.convolve_mirror <- function(x, kernel) {
  n <- length(x)
  half_k <- (length(kernel) - 1L) %/% 2L

  if (half_k == 0L || n <= 1L) return(x)

  # Mirror-pad edges
  pad_left <- rev(x[seq_len(min(half_k, n))])
  if (length(pad_left) < half_k) {
    pad_left <- rep_len(pad_left, half_k)
  }
  pad_right <- rev(x[seq(max(1L, n - half_k + 1L), n)])
  if (length(pad_right) < half_k) {
    pad_right <- rep_len(pad_right, half_k)
  }
  padded <- c(pad_left, x, pad_right)

  # Manual convolution (handles complex values)
  n_padded <- length(padded)
  result <- vector(mode = if (is.complex(x)) "complex" else "numeric",
                   length = n)
  for (i in seq_len(n)) {
    idx <- (i - 1L) + seq_along(kernel)
    result[i] <- sum(padded[idx] * kernel)
  }
  result
}
