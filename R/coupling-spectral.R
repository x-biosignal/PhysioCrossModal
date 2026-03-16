# Spectral coupling functions for cross-modal analysis
#
# Functions for computing spectral coherence and cross-spectral density
# between signals from different physiological modalities.

#' Magnitude-squared coherence between two signals
#'
#' Computes the magnitude-squared coherence between two signals, which
#' quantifies the linear frequency-domain relationship between them.
#' The coherence is defined as:
#' \deqn{C_{xy}(f) = \frac{|S_{xy}(f)|^2}{S_{xx}(f) \cdot S_{yy}(f)}}
#' where \eqn{S_{xy}} is the cross-spectral density and \eqn{S_{xx}},
#' \eqn{S_{yy}} are the auto-spectral densities.
#'
#' Unlike the `coherence()` function in PhysioAnalysis (which computes
#' coherence across all channel pairs within a single PhysioExperiment),
#' this function is designed for cross-modal analysis: it takes two separate
#' signals (or PhysioExperiment / MultiPhysioExperiment objects) and computes
#' the pairwise coherence between them.
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when `x` is an MPE.
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param freq_range Numeric vector of length 2 giving the frequency range
#'   \code{c(low, high)} in Hz to restrict the output. NULL returns all
#'   frequencies.
#' @param nperseg Integer segment length for Welch's method (default 256).
#' @param noverlap Integer overlap between segments, or NULL for
#'   \code{floor(nperseg / 2)}.
#' @param method Character spectral estimation method. Currently only
#'   \code{"welch"} is implemented; \code{"multitaper"} is reserved for
#'   future use.
#' @param modality_x Character modality name in MPE for the x signal.
#' @param modality_y Character modality name in MPE for the y signal.
#' @param channels_x Integer which channel to extract from x (default 1).
#' @param channels_y Integer which channel to extract from y (default 1).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{coherence}{Numeric vector of magnitude-squared coherence values
#'       in \eqn{[0, 1]} at each frequency.}
#'     \item{frequencies}{Numeric vector of corresponding frequencies in Hz.}
#'     \item{confidence_limit}{Numeric scalar giving the 95\% confidence
#'       threshold: \code{1 - alpha^(1 / (L - 1))} where alpha = 0.05 and
#'       L is the number of segments (Halliday et al. 1995). Coherence
#'       values above this limit are statistically significant.}
#'     \item{n_segments}{Integer number of segments used in Welch estimation.}
#'   }
#'
#' @references
#' Carter, G. C. (1987). Coherence and time delay estimation.
#' \emph{Proceedings of the IEEE}, 75(2), 236--255.
#'
#' Halliday, D. M., Rosenberg, J. R., Amjad, A. M., Breeze, P.,
#' Conway, B. A., & Farmer, S. F. (1995). A framework for the analysis of
#' mixed time series/point process data -- theory and application to the
#' study of physiological tremor, single motor unit discharges and
#' electromyograms. \emph{Progress in Biophysics and Molecular Biology},
#' 64(2--3), 237--278.
#'
#' @seealso [crossSpectrum()], [multitaperCoherence()],
#'   [coherenceMatrix()], [couplingAnalysis()]
#' @export
#' @examples
#' # Two coupled 20 Hz sinusoids
#' sr <- 500
#' t <- seq(0, 10, length.out = sr * 10)
#' x <- sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
#' y <- 0.8 * sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
#' result <- coherence(x, y, sr = sr)
#' plot(result$frequencies, result$coherence, type = "l",
#'      xlab = "Frequency (Hz)", ylab = "Coherence")
#' abline(h = result$confidence_limit, lty = 2, col = "red")
coherence <- function(x, y = NULL, sr = NULL, freq_range = NULL,
                      nperseg = 256L, noverlap = NULL,
                      method = c("welch", "multitaper"),
                      modality_x = NULL, modality_y = NULL,
                      channels_x = 1L, channels_y = 1L) {

  method <- match.arg(method)
  nperseg <- as.integer(nperseg)

  if (method == "multitaper") {
    stop("method = 'multitaper' is not yet implemented; use 'welch'",
         call. = FALSE)
  }

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

  # Compute auto- and cross-spectral densities
  psd_xx <- .welch_psd(sig_x, nperseg = nperseg, noverlap = noverlap,
                        sr = sr_val)
  psd_yy <- .welch_psd(sig_y, nperseg = nperseg, noverlap = noverlap,
                        sr = sr_val)
  csd_xy <- .welch_csd(sig_x, sig_y, nperseg = nperseg, noverlap = noverlap,
                        sr = sr_val)

  # Both .welch_psd() and .welch_csd() now apply one-sided spectrum
  # doubling consistently, so the 2x factors cancel in the coherence ratio:
  #   |2*Sxy|^2 / (2*Sxx * 2*Syy) = |Sxy|^2 / (Sxx*Syy)
  csd_vals <- csd_xy$csd

  # Magnitude-squared coherence: |Sxy(f)|^2 / (Sxx(f) * Syy(f))
  coh <- (Mod(csd_vals)^2) / (psd_xx$psd * psd_yy$psd)

  # Clamp to [0, 1] to handle numerical imprecision
  coh <- pmin(pmax(coh, 0), 1)

  freqs <- csd_xy$frequencies

  # Compute number of segments used (same logic as .welch_psd)
  n <- min_len
  effective_nperseg <- nperseg
  effective_noverlap <- if (is.null(noverlap)) floor(nperseg / 2L) else as.integer(noverlap)
  if (n < effective_nperseg) {
    effective_nperseg <- n
    effective_noverlap <- 0L
  }
  step <- effective_nperseg - effective_noverlap
  n_segments <- 0L
  for (seg in seq_len(max(1L, floor((n - effective_noverlap) / step)))) {
    start <- (seg - 1L) * step + 1L
    end_pos <- start + effective_nperseg - 1L
    if (end_pos > n) break
    n_segments <- n_segments + 1L
  }

  # Confidence limit: 1 - alpha^(1 / (L - 1))
  # where alpha = 0.05 and L = number of segments (Halliday et al. 1995)
  if (n_segments > 1L) {
    confidence_limit <- 1 - 0.05^(1 / (n_segments - 1))
  } else {
    confidence_limit <- 1.0  # Cannot determine significance with 1 segment
  }

  # Apply frequency range filter if specified
  if (!is.null(freq_range)) {
    stopifnot(is.numeric(freq_range), length(freq_range) == 2)
    freq_idx <- which(freqs >= freq_range[1] & freqs <= freq_range[2])
    freqs <- freqs[freq_idx]
    coh <- coh[freq_idx]
  }

  list(
    coherence = coh,
    frequencies = freqs,
    confidence_limit = confidence_limit,
    n_segments = n_segments
  )
}


#' Cross-spectral density between two signals
#'
#' Computes the cross-spectral density (CSD) between two signals using
#' Welch's method. The CSD captures both the magnitude and phase relationship
#' between two signals as a function of frequency.
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when `x` is an MPE.
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param nperseg Integer segment length for Welch's method (default 256).
#' @param noverlap Integer overlap between segments, or NULL for
#'   \code{floor(nperseg / 2)}.
#' @param modality_x Character modality name in MPE for the x signal.
#' @param modality_y Character modality name in MPE for the y signal.
#' @param channels_x Integer which channel to extract from x (default 1).
#' @param channels_y Integer which channel to extract from y (default 1).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{csd}{Complex vector of cross-spectral density values.}
#'     \item{frequencies}{Numeric vector of corresponding frequencies in Hz.}
#'     \item{magnitude}{Numeric vector of CSD magnitude (\code{Mod(csd)}).}
#'     \item{phase}{Numeric vector of CSD phase in radians (\code{Arg(csd)}).}
#'   }
#'
#' @references
#' Carter, G. C. (1987). Coherence and time delay estimation.
#' \emph{Proceedings of the IEEE}, 75(2), 236--255.
#'
#' Halliday, D. M., Rosenberg, J. R., Amjad, A. M., Breeze, P.,
#' Conway, B. A., & Farmer, S. F. (1995). A framework for the analysis of
#' mixed time series/point process data -- theory and application to the
#' study of physiological tremor, single motor unit discharges and
#' electromyograms. \emph{Progress in Biophysics and Molecular Biology},
#' 64(2--3), 237--278.
#'
#' @seealso [coherence()], [multitaperCoherence()], [couplingAnalysis()]
#' @export
#' @examples
#' sr <- 500
#' t <- seq(0, 5, length.out = sr * 5)
#' x <- sin(2 * pi * 10 * t) + 0.1 * rnorm(length(t))
#' y <- cos(2 * pi * 10 * t) + 0.1 * rnorm(length(t))
#' result <- crossSpectrum(x, y, sr = sr)
#' plot(result$frequencies, result$magnitude, type = "l",
#'      xlab = "Frequency (Hz)", ylab = "|CSD|")
crossSpectrum <- function(x, y = NULL, sr = NULL,
                          nperseg = 256L, noverlap = NULL,
                          modality_x = NULL, modality_y = NULL,
                          channels_x = 1L, channels_y = 1L) {

  nperseg <- as.integer(nperseg)

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

  # Compute cross-spectral density
  result <- .welch_csd(sig_x, sig_y, nperseg = nperseg, noverlap = noverlap,
                        sr = sr_val)

  list(
    csd = result$csd,
    frequencies = result$frequencies,
    magnitude = Mod(result$csd),
    phase = Arg(result$csd)
  )
}


#' Multitaper coherence
#'
#' Computes magnitude-squared coherence using multitaper spectral estimation
#' with Discrete Prolate Spheroidal Sequences (DPSS / Slepian tapers).
#' Multitaper estimation provides lower variance than Welch's method for
#' short signals or when frequency resolution is critical.
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when \code{x} is an
#'   MPE.
#' @param sr Numeric sampling rate in Hz (required when x/y are numeric).
#' @param nw Numeric time-bandwidth product (default 4). Controls the
#'   trade-off between frequency resolution and variance.
#' @param k Integer number of tapers (default \code{2 * nw - 1}).
#' @param freq_range Optional numeric vector \code{c(low, high)} to restrict
#'   output frequencies.
#' @param modality_x,modality_y Character modality names for MPE input.
#' @param channels_x,channels_y Integer channel indices (default 1).
#' @param ... Currently unused.
#'
#' @return A list with components:
#'   \describe{
#'     \item{coherence}{Numeric vector of magnitude-squared coherence values
#'       in \eqn{[0, 1]}.}
#'     \item{frequencies}{Numeric vector of corresponding frequencies in Hz.}
#'     \item{confidence_limit}{Numeric scalar giving the 95\% confidence
#'       threshold based on the number of tapers.}
#'   }
#'
#' @references
#' Thomson, D. J. (1982). Spectrum estimation and harmonic analysis.
#' \emph{Proceedings of the IEEE}, 70(9), 1055--1096.
#'
#' Carter, G. C. (1987). Coherence and time delay estimation.
#' \emph{Proceedings of the IEEE}, 75(2), 236--255.
#'
#' @seealso [coherence()], [crossSpectrum()], [couplingAnalysis()]
#' @export
#' @examples
#' sr <- 500
#' t <- seq(0, 10, length.out = sr * 10)
#' x <- sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
#' y <- 0.8 * sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
#' result <- multitaperCoherence(x, y, sr = sr)
multitaperCoherence <- function(x, y = NULL, sr = NULL,
                                 nw = 4, k = NULL,
                                 freq_range = NULL,
                                 modality_x = NULL, modality_y = NULL,
                                 channels_x = 1L, channels_y = 1L, ...) {

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

  # Compute multitaper CSD
  mt <- .multitaper_csd(sig_x, sig_y, sr_val, nw = nw, k = k)

  # Magnitude-squared coherence: |Sxy|^2 / (Sxx * Syy)
  denom <- mt$psd_x * mt$psd_y
  coh <- ifelse(denom > .Machine$double.eps,
                Mod(mt$csd)^2 / denom,
                0)
  coh <- pmin(pmax(coh, 0), 1)

  freqs <- mt$frequencies

  # Confidence limit based on number of tapers
  if (is.null(k)) k <- as.integer(2 * nw - 1)
  k <- min(k, n)
  if (k > 1L) {
    confidence_limit <- 1 - 0.05^(1 / (k - 1))
  } else {
    confidence_limit <- 1.0
  }

  # Apply frequency range filter
  if (!is.null(freq_range)) {
    stopifnot(is.numeric(freq_range), length(freq_range) == 2)
    freq_idx <- which(freqs >= freq_range[1] & freqs <= freq_range[2])
    freqs <- freqs[freq_idx]
    coh <- coh[freq_idx]
  }

  list(
    coherence = coh,
    frequencies = freqs,
    confidence_limit = confidence_limit
  )
}
