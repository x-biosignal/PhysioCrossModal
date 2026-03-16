# Phase coupling functions for cross-modal analysis
#
# Implements Phase Locking Value (PLV), Phase Lag Index (PLI), and
# weighted Phase Lag Index (wPLI) for quantifying phase synchrony
# between physiological signals from different modalities.

# ---- Phase Locking Value (PLV) -----------------------------------------------

#' Phase Locking Value (PLV)
#'
#' Computes the Phase Locking Value between two signals, measuring the
#' consistency of phase difference across time. A PLV of 1 indicates
#' perfect phase locking; a PLV of 0 indicates no consistent phase
#' relationship.
#'
#' The signals are bandpass-filtered to `freq_band` and then Hilbert-
#' transformed to extract instantaneous phase. PLV is computed as:
#' \deqn{\text{PLV} = \left|\frac{1}{N} \sum_{t=1}^{N} e^{i(\phi_x(t) - \phi_y(t))}\right|}
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment (NULL when `x` is an MPE).
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param freq_band Numeric vector of length 2 specifying the frequency band
#'   \code{c(low, high)} in Hz. Required.
#' @param method Character; phase extraction method. Currently only
#'   `"hilbert"` is supported. `"wavelet"` is reserved for future use.
#' @param modality_x Character modality name in MPE for x signal.
#' @param modality_y Character modality name in MPE for y signal.
#' @param channels_x Integer channel index to extract from x (default 1).
#' @param channels_y Integer channel index to extract from y (default 1).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{plv}{Numeric scalar in \[0, 1\] -- the Phase Locking Value.}
#'     \item{phase_diff}{Numeric vector of instantaneous phase differences
#'       (radians, \[-pi, pi\]).}
#'   }
#'
#' @references
#' Lachaux, J.-P., Rodriguez, E., Martinerie, J., & Varela, F. J. (1999).
#' Measuring phase synchrony in brain signals. \emph{Human Brain Mapping},
#' 8(4), 194--208.
#'
#' @seealso [phaseLagIndex()], [weightedPLI()], [waveletPLV()],
#'   [couplingAnalysis()]
#' @export
#' @examples
#' # Two phase-locked sinusoids
#' sr <- 500
#' t <- seq(0, 2, length.out = sr * 2)
#' x <- sin(2 * pi * 10 * t)
#' y <- sin(2 * pi * 10 * t + pi / 4)
#' result <- phaseLockingValue(x, y, sr = sr, freq_band = c(8, 12))
#' result$plv
phaseLockingValue <- function(x, y = NULL, sr = NULL, freq_band,
                              method = c("hilbert", "wavelet"),
                              modality_x = NULL, modality_y = NULL,
                              channels_x = 1L, channels_y = 1L) {

  # Validate freq_band
  if (missing(freq_band)) {
    stop("'freq_band' is required: a numeric vector of length 2 (c(low, high))",
         call. = FALSE)
  }
  if (!is.numeric(freq_band) || length(freq_band) != 2) {
    stop("'freq_band' must be a numeric vector of length 2", call. = FALSE)
  }

  method <- match.arg(method)
  if (method == "wavelet") {
    stop("method = 'wavelet' is not yet implemented", call. = FALSE)
  }

  # Extract signal pair
  pair <- .extract_signal_pair(x, y, sr = sr,
                               modality_x = modality_x,
                               modality_y = modality_y,
                               channels_x = channels_x,
                               channels_y = channels_y)
  sig_x <- pair$x
  sig_y <- pair$y
  sr_use <- pair$sr

  # Bandpass filter both signals to freq_band
  filt_x <- .bandpass_filter(sig_x, low = freq_band[1], high = freq_band[2],
                             sr = sr_use)
  filt_y <- .bandpass_filter(sig_y, low = freq_band[1], high = freq_band[2],
                             sr = sr_use)

  # Extract instantaneous phase via Hilbert transform
  phase_x <- .hilbert_phase(filt_x)
  phase_y <- .hilbert_phase(filt_y)

  # Compute phase difference
  phase_diff <- phase_x - phase_y
  # Wrap to [-pi, pi]
  phase_diff <- atan2(sin(phase_diff), cos(phase_diff))

  # PLV = |mean(exp(i * phase_diff))|
  plv <- abs(mean(exp(complex(imaginary = phase_diff))))

  list(plv = plv, phase_diff = phase_diff)
}

# ---- Phase Lag Index (PLI) ---------------------------------------------------

#' Phase Lag Index (PLI)
#'
#' Computes the Phase Lag Index between two signals. PLI measures the
#' asymmetry of the distribution of phase differences and is insensitive
#' to zero-lag (volume conduction) effects.
#'
#' \deqn{\text{PLI} = \left|\frac{1}{N} \sum_{t=1}^{N} \text{sign}(\sin(\phi_x(t) - \phi_y(t)))\right|}
#'
#' A PLI of 0 indicates either no coupling or symmetric (zero-lag) coupling.
#' A PLI of 1 indicates perfect non-zero-lag phase locking.
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment (NULL when `x` is an MPE).
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param freq_band Numeric vector of length 2 specifying the frequency band
#'   \code{c(low, high)} in Hz. Required.
#' @param modality_x Character modality name in MPE for x signal.
#' @param modality_y Character modality name in MPE for y signal.
#' @param channels_x Integer channel index to extract from x (default 1).
#' @param channels_y Integer channel index to extract from y (default 1).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{pli}{Numeric scalar in \[0, 1\] -- the Phase Lag Index.}
#'     \item{phase_diff}{Numeric vector of instantaneous phase differences
#'       (radians).}
#'   }
#'
#' @references
#' Stam, C. J., Nolte, G., & Daffertshofer, A. (2007). Phase lag index:
#' Assessment of functional connectivity from multi channel EEG and MEG
#' with diminished bias from common sources. \emph{Human Brain Mapping},
#' 28(11), 1178--1193.
#'
#' @seealso [phaseLockingValue()], [weightedPLI()], [couplingAnalysis()]
#' @export
#' @examples
#' sr <- 500
#' t <- seq(0, 2, length.out = sr * 2)
#' x <- sin(2 * pi * 10 * t)
#' y <- sin(2 * pi * 10 * t + pi / 3)
#' result <- phaseLagIndex(x, y, sr = sr, freq_band = c(8, 12))
#' result$pli
phaseLagIndex <- function(x, y = NULL, sr = NULL, freq_band,
                          modality_x = NULL, modality_y = NULL,
                          channels_x = 1L, channels_y = 1L) {

  # Validate freq_band
  if (missing(freq_band)) {
    stop("'freq_band' is required: a numeric vector of length 2 (c(low, high))",
         call. = FALSE)
  }
  if (!is.numeric(freq_band) || length(freq_band) != 2) {
    stop("'freq_band' must be a numeric vector of length 2", call. = FALSE)
  }

  # Extract signal pair
  pair <- .extract_signal_pair(x, y, sr = sr,
                               modality_x = modality_x,
                               modality_y = modality_y,
                               channels_x = channels_x,
                               channels_y = channels_y)
  sig_x <- pair$x
  sig_y <- pair$y
  sr_use <- pair$sr

  # Bandpass filter both signals to freq_band
  filt_x <- .bandpass_filter(sig_x, low = freq_band[1], high = freq_band[2],
                             sr = sr_use)
  filt_y <- .bandpass_filter(sig_y, low = freq_band[1], high = freq_band[2],
                             sr = sr_use)

  # Extract instantaneous phase via Hilbert transform
  phase_x <- .hilbert_phase(filt_x)
  phase_y <- .hilbert_phase(filt_y)

  # Compute phase difference
  phase_diff <- phase_x - phase_y

  # PLI = |mean(sign(sin(phase_diff)))|
  pli <- abs(mean(sign(sin(phase_diff))))

  list(pli = pli, phase_diff = phase_diff)
}

# ---- Weighted Phase Lag Index (wPLI) -----------------------------------------

#' Weighted Phase Lag Index (wPLI)
#'
#' Computes the weighted Phase Lag Index between two signals. wPLI weights
#' each phase-difference sample by the magnitude of the imaginary component
#' of the cross-spectrum, reducing the influence of noise sources that produce
#' phase differences near 0 or pi.
#'
#' \deqn{\text{wPLI} = \frac{\left|\sum |\text{Im}(S_{xy})| \cdot \text{sign}(\text{Im}(S_{xy}))\right|}{\sum |\text{Im}(S_{xy})|}}{wPLI = |sum(|Im(Sxy)| * sign(Im(Sxy)))| / sum(|Im(Sxy)|)}
#'
#' When `debiased = TRUE`, the debiased estimator from Vinck et al. (2011) is
#' also returned:
#' \deqn{\text{wPLI}_{\text{debiased}} = \frac{(\sum \text{Im}(S_{xy}))^2 - \sum \text{Im}(S_{xy})^2}{(\sum |\text{Im}(S_{xy})|)^2 - \sum \text{Im}(S_{xy})^2}}
#'
#' @references
#' Vinck, M., Oostenveld, R., van Wingerden, M., Battaglia, F., &
#' Pennartz, C. M. A. (2011). An improved index of phase-synchronization
#' for electrophysiological data in the presence of volume-conduction, noise
#' and sample-size bias. \emph{NeuroImage}, 55(4), 1548--1565.
#'
#' Lachaux, J.-P., Rodriguez, E., Martinerie, J., & Varela, F. J. (1999).
#' Measuring phase synchrony in brain signals. \emph{Human Brain Mapping},
#' 8(4), 194--208.
#'
#' @seealso [phaseLockingValue()], [phaseLagIndex()], [couplingAnalysis()]
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment (NULL when `x` is an MPE).
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param freq_band Numeric vector of length 2 specifying the frequency band
#'   \code{c(low, high)} in Hz. Required.
#' @param debiased Logical; if TRUE (default), also compute the debiased wPLI
#'   estimator.
#' @param modality_x Character modality name in MPE for x signal.
#' @param modality_y Character modality name in MPE for y signal.
#' @param channels_x Integer channel index to extract from x (default 1).
#' @param channels_y Integer channel index to extract from y (default 1).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{wpli}{Numeric scalar in \[0, 1\] -- the weighted PLI.}
#'     \item{wpli_debiased}{Numeric scalar -- the debiased wPLI, or NULL if
#'       `debiased = FALSE`.}
#'     \item{n_samples}{Integer -- number of time-domain samples used.}
#'   }
#'
#' @export
#' @examples
#' sr <- 500
#' t <- seq(0, 2, length.out = sr * 2)
#' x <- sin(2 * pi * 10 * t)
#' y <- sin(2 * pi * 10 * t + pi / 4)
#' result <- weightedPLI(x, y, sr = sr, freq_band = c(8, 12))
#' result$wpli
weightedPLI <- function(x, y = NULL, sr = NULL, freq_band, debiased = TRUE,
                        modality_x = NULL, modality_y = NULL,
                        channels_x = 1L, channels_y = 1L) {

  # Validate freq_band
  if (missing(freq_band)) {
    stop("'freq_band' is required: a numeric vector of length 2 (c(low, high))",
         call. = FALSE)
  }
  if (!is.numeric(freq_band) || length(freq_band) != 2) {
    stop("'freq_band' must be a numeric vector of length 2", call. = FALSE)
  }

  # Extract signal pair
  pair <- .extract_signal_pair(x, y, sr = sr,
                               modality_x = modality_x,
                               modality_y = modality_y,
                               channels_x = channels_x,
                               channels_y = channels_y)
  sig_x <- pair$x
  sig_y <- pair$y
  sr_use <- pair$sr

  # Bandpass filter both signals to freq_band
  filt_x <- .bandpass_filter(sig_x, low = freq_band[1], high = freq_band[2],
                             sr = sr_use)
  filt_y <- .bandpass_filter(sig_y, low = freq_band[1], high = freq_band[2],
                             sr = sr_use)

  # Compute analytic signals via Hilbert transform
  analytic_x <- .hilbert_analytic(filt_x)
  analytic_y <- .hilbert_analytic(filt_y)

  # Cross-spectrum at each time point
  cross <- analytic_x * Conj(analytic_y)
  imag_cross <- Im(cross)

  n_samples <- length(imag_cross)

  # wPLI = |sum(|Im(Sxy)| * sign(Im(Sxy)))| / sum(|Im(Sxy)|)
  abs_imag <- abs(imag_cross)
  sum_abs_imag <- sum(abs_imag)

  if (sum_abs_imag < .Machine$double.eps) {
    # Degenerate case: no imaginary component (e.g. zero-lag or DC)
    wpli_val <- 0
  } else {
    wpli_val <- abs(sum(abs_imag * sign(imag_cross))) / sum_abs_imag
  }

  # Debiased wPLI (Vinck et al. 2011, NeuroImage 55:1548-1565, Eq. 6):
  #   (sum(Im(Sxy))^2 - sum(Im(Sxy)^2)) / (sum(|Im(Sxy)|)^2 - sum(Im(Sxy)^2))
  wpli_debiased <- NULL
  if (debiased) {
    if (n_samples > 1) {
      sum_imag <- sum(imag_cross)
      sum_imag_sq <- sum(imag_cross^2)
      sum_abs_sq <- sum_abs_imag^2
      denom <- sum_abs_sq - sum_imag_sq
      if (abs(denom) < .Machine$double.eps) {
        wpli_debiased <- 0
      } else {
        wpli_debiased <- (sum_imag^2 - sum_imag_sq) / denom
      }
    } else {
      wpli_debiased <- NA_real_
    }
  }

  list(
    wpli = wpli_val,
    wpli_debiased = wpli_debiased,
    n_samples = as.integer(n_samples)
  )
}
