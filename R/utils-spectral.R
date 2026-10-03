# Shared spectral utilities for PhysioCrossModal
#
# Internal (non-exported) helper functions used by all coupling modules:
# - Welch PSD / CSD estimation
# - Hilbert transform (analytic signal, phase, envelope)
# - Bandpass filtering (Butterworth + FFT fallback)
# - Window functions
# - Signal-pair extraction from various input types

# ---- Window Functions -------------------------------------------------------

#' Hanning (raised cosine) window
#'
#' Generates a symmetric Hanning window of length `n`:
#' \code{0.5 * (1 - cos(2 * pi * k / (n - 1)))}.
#'
#' @param n Integer window length.
#' @return Numeric vector of length `n`.
#' @noRd
.window_hanning <- function(n) {

  n <- as.integer(n)
  if (n <= 0L) return(numeric(0L))
  if (n == 1L) return(1)
  0.5 * (1 - cos(2 * pi * seq(0L, n - 1L) / (n - 1L)))
}

#' Hamming window
#'
#' Generates a Hamming window of length `n`:
#' \code{0.54 - 0.46 * cos(2 * pi * k / (n - 1))}.
#'
#' @param n Integer window length.
#' @return Numeric vector of length `n`.
#' @noRd
.window_hamming <- function(n) {
  n <- as.integer(n)
  if (n <= 0L) return(numeric(0L))
  if (n == 1L) return(1)
  0.54 - 0.46 * cos(2 * pi * seq(0L, n - 1L) / (n - 1L))
}

#' Resolve a window function by name
#'
#' @param window Character name or numeric vector.
#' @param n Integer window length (ignored when `window` is already numeric).
#' @return Numeric vector of length `n`.
#' @noRd
.resolve_window <- function(window, n) {
  if (is.numeric(window)) {
    if (length(window) != n) {
      stop("Numeric window must have length equal to nperseg (", n, ")")
    }
    return(window)
  }
  window <- tolower(window)
  switch(window,
    hanning = .window_hanning(n),
    hann    = .window_hanning(n),
    hamming = .window_hamming(n),
    stop("Unknown window type: '", window, "'. Use 'hanning' or 'hamming'.")
  )
}

# ---- Welch PSD Estimation ---------------------------------------------------

#' Welch PSD estimation
#'
#' Estimates the one-sided power spectral density of a real-valued signal using
#' overlapping windowed segments (Welch's method).
#'
#' @param x Numeric vector — the input signal.
#' @param nperseg Integer — length of each FFT segment (default 256).
#' @param noverlap Integer or NULL — overlap between segments. Defaults to
#'   `floor(nperseg / 2)`.
#' @param sr Numeric — sampling rate in Hz (default 1).
#' @param window Character name of window function or numeric vector.
#' @return Named list with components:
#'   \describe{
#'     \item{psd}{Numeric vector of PSD values (one-sided).}
#'     \item{frequencies}{Numeric vector of corresponding frequencies in Hz.}
#'   }
#' @noRd
.welch_psd <- function(x, nperseg = 256L, noverlap = NULL, sr = 1,
                       window = "hanning") {
  x <- as.numeric(x)
  n <- length(x)
  nperseg <- as.integer(nperseg)
  if (is.null(noverlap)) noverlap <- floor(nperseg / 2L)
  noverlap <- as.integer(noverlap)

  # Validation
  if (nperseg < 2L) stop("nperseg must be >= 2")
  if (noverlap < 0L || noverlap >= nperseg) {
    stop("noverlap must be in [0, nperseg - 1]")
  }

  # If signal is shorter than one segment, adapt
  if (n < nperseg) {
    nperseg <- n
    noverlap <- 0L
  }

  win <- .resolve_window(window, nperseg)
  step <- nperseg - noverlap
  n_segments <- max(1L, floor((n - noverlap) / step))
  n_freqs <- floor(nperseg / 2) + 1L

  psd <- numeric(n_freqs)
  actual_segments <- 0L


  for (seg in seq_len(n_segments)) {
    start <- (seg - 1L) * step + 1L
    end <- start + nperseg - 1L
    if (end > n) break

    segment <- x[start:end] * win
    fft_result <- stats::fft(segment)
    psd <- psd + Mod(fft_result[seq_len(n_freqs)])^2
    actual_segments <- actual_segments + 1L
  }

  if (actual_segments == 0L) {
    stop("No complete segments could be formed from the signal")
  }

  # Normalize: divide by number of segments and by (sr * sum(window^2))
  psd <- psd / actual_segments
  psd <- psd / (sr * sum(win^2))

  # Double the non-DC, non-Nyquist bins for one-sided spectrum
  if (n_freqs > 2L) {
    psd[2:(n_freqs - 1L)] <- 2 * psd[2:(n_freqs - 1L)]
  }

  list(
    psd = psd,
    frequencies = seq(0, sr / 2, length.out = n_freqs)
  )
}

# ---- Cross-Spectral Density ------------------------------------------------

#' Cross-spectral density via Welch method
#'
#' Computes the one-sided cross-spectral density between two equal-length
#' signals using the Welch approach.
#'
#' @param x Numeric vector — first signal.
#' @param y Numeric vector — second signal (same length as `x`).
#' @param nperseg Integer — segment length (default 256).
#' @param noverlap Integer or NULL — overlap (default `nperseg / 2`).
#' @param sr Numeric — sampling rate in Hz (default 1).
#' @param window Character name of window or numeric vector.
#' @return Named list with components:
#'   \describe{
#'     \item{csd}{Complex vector — cross-spectral density.}
#'     \item{frequencies}{Numeric vector of corresponding frequencies.}
#'   }
#' @noRd
.welch_csd <- function(x, y, nperseg = 256L, noverlap = NULL, sr = 1,
                       window = "hanning") {
  x <- as.numeric(x)
  y <- as.numeric(y)
  n <- length(x)
  if (length(y) != n) {
    stop("x and y must have the same length")
  }

  nperseg <- as.integer(nperseg)
  if (is.null(noverlap)) noverlap <- floor(nperseg / 2L)
  noverlap <- as.integer(noverlap)

  if (nperseg < 2L) stop("nperseg must be >= 2")
  if (noverlap < 0L || noverlap >= nperseg) {
    stop("noverlap must be in [0, nperseg - 1]")
  }

  if (n < nperseg) {
    nperseg <- n
    noverlap <- 0L
  }

  win <- .resolve_window(window, nperseg)
  step <- nperseg - noverlap
  n_segments <- max(1L, floor((n - noverlap) / step))
  n_freqs <- floor(nperseg / 2) + 1L

  csd <- complex(n_freqs)
  actual_segments <- 0L

  for (seg in seq_len(n_segments)) {
    start <- (seg - 1L) * step + 1L
    end <- start + nperseg - 1L
    if (end > n) break

    seg_x <- x[start:end] * win
    seg_y <- y[start:end] * win

    fft_x <- stats::fft(seg_x)
    fft_y <- stats::fft(seg_y)

    # Cross-spectrum: X * conj(Y)
    csd <- csd + fft_x[seq_len(n_freqs)] * Conj(fft_y[seq_len(n_freqs)])
    actual_segments <- actual_segments + 1L
  }

  if (actual_segments == 0L) {
    stop("No complete segments could be formed from the signals")
  }

  csd <- csd / actual_segments
  csd <- csd / (sr * sum(win^2))

  # Apply one-sided doubling to match .welch_psd() convention
  if (n_freqs > 2L) {
    csd[2:(n_freqs - 1L)] <- 2 * csd[2:(n_freqs - 1L)]
  }

  list(
    csd = csd,
    frequencies = seq(0, sr / 2, length.out = n_freqs)
  )
}

# ---- Hilbert Transform -----------------------------------------------------

#' Hilbert transform — analytic signal
#'
#' Computes the analytic signal via FFT-based Hilbert transform. The result is
#' a complex vector whose real part equals the original signal and whose
#' imaginary part is the Hilbert transform.
#'
#' @param x Numeric vector.
#' @return Complex vector of same length as `x`.
#' @noRd
.hilbert_analytic <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0L) return(complex(0L))

  fft_x <- stats::fft(x)

  # Build the analytic-signal multiplier h:
  #   h[1] = 1 (DC), h[N/2+1] = 1 (Nyquist for even N),
  #   h[2..N/2] = 2 (positive frequencies), rest = 0
  h <- numeric(n)
  h[1] <- 1
  if (n %% 2 == 0) {
    if (n >= 4L) h[2:(n / 2)] <- 2
    h[n / 2 + 1] <- 1  # Nyquist
  } else {
    if (n >= 3L) h[2:((n + 1) / 2)] <- 2
  }

  stats::fft(fft_x * h, inverse = TRUE) / n
}

#' Instantaneous phase via Hilbert transform
#'
#' Returns the instantaneous phase (in radians, range \eqn{[-\pi, \pi]}) of the
#' analytic signal computed from `x`.
#'
#' @param x Numeric vector.
#' @return Numeric vector of phase angles (radians).
#' @noRd
.hilbert_phase <- function(x) {
  analytic <- .hilbert_analytic(x)
  Arg(analytic)
}

#' Instantaneous amplitude envelope via Hilbert transform
#'
#' Returns the envelope (instantaneous amplitude) of the signal, i.e.
#' `Mod(analytic_signal)`.
#'
#' @param x Numeric vector.
#' @return Numeric vector of amplitude values.
#' @noRd
.hilbert_envelope <- function(x) {
  analytic <- .hilbert_analytic(x)
  Mod(analytic)
}

# ---- Bandpass Filter --------------------------------------------------------

#' Bandpass filter (Butterworth or FIR, with FFT fallback)
#'
#' Applies a bandpass filter to a signal. First tries to use the `signal`
#' package for an IIR (Butterworth) or FIR filter with zero-phase filtering
#' (filtfilt). Falls back to a simple FFT-based brick-wall filter when the
#' `signal` package is not available.
#'
#' @param x Numeric vector — input signal.
#' @param low Numeric — low cutoff frequency in Hz.
#' @param high Numeric — high cutoff frequency in Hz.
#' @param sr Numeric — sampling rate in Hz.
#' @param order Integer — filter order (default 4).
#' @param type Character — filter type: `"butter"` (IIR Butterworth, default)
#'   or `"fir"` (FIR via `signal::fir1`).
#' @param clamp Logical; clamp normalized cutoffs to `[0.01, 0.99]`. Existing
#'   callers retain this numerical safeguard. Audited callers may set it to
#'   `FALSE` after validating exact pass bands.
#' @return Numeric vector — filtered signal (same length as `x`).
#' @noRd
.bandpass_filter <- function(x, low, high, sr, order = 4L,
                             type = c("butter", "fir"), clamp = TRUE) {
  type <- match.arg(type)
  x <- as.numeric(x)
  n <- length(x)
  nyq <- sr / 2

  # Input validation
  if (low <= 0 || high <= 0) stop("Cutoff frequencies must be positive")
  if (low >= high) stop("'low' must be less than 'high'")
  if (high >= nyq) {
    stop(sprintf(
      "high cutoff (%.1f Hz) must be below Nyquist frequency (%.1f Hz)",
      high, nyq
    ))
  }

  Wn <- c(low / nyq, high / nyq)
  if (isTRUE(clamp)) {
    Wn_clamped <- pmax(pmin(Wn, 0.99), 0.01)
    if (!isTRUE(all.equal(Wn, Wn_clamped))) {
      warning(sprintf(
        "Normalized cutoff frequencies clamped from [%.4f, %.4f] to [%.4f, %.4f] for numerical stability",
        Wn[1], Wn[2], Wn_clamped[1], Wn_clamped[2]
      ), call. = FALSE)
    }
    Wn <- Wn_clamped
  }

  tryCatch({
    if (!requireNamespace("signal", quietly = TRUE)) {
      stop("signal package not available")
    }
    if (type == "butter") {
      bf <- signal::butter(order, Wn, type = "pass")
      signal::filtfilt(bf, x)
    } else {
      # FIR filter via fir1
      fir_order <- 2L * order * 10L  # reasonable FIR order
      bf <- signal::fir1(fir_order, Wn, type = "pass")
      signal::filtfilt(bf, x)
    }
  }, error = function(e) {
    # Fallback: FFT-based brick-wall bandpass filter
    fft_x <- stats::fft(x)
    freqs <- seq(0, sr - sr / n, length.out = n)

    # Create symmetric bandpass mask (positive + negative freqs)
    mask <- (freqs >= low & freqs <= high) |
            (freqs >= (sr - high) & freqs <= (sr - low))
    fft_x[!mask] <- 0

    Re(stats::fft(fft_x, inverse = TRUE)) / n
  })
}

# ---- Low-pass Filter (Anti-aliasing) ----------------------------------------

#' Low-pass filter for anti-aliasing before downsampling
#'
#' Applies a low-pass filter at the specified cutoff frequency. Uses the
#' `signal` package when available; falls back to FFT brick-wall filtering.
#'
#' @param x Numeric vector — input signal.
#' @param cutoff Numeric — cutoff frequency in Hz.
#' @param sr Numeric — sampling rate in Hz.
#' @param order Integer — filter order (default 8).
#' @return Numeric vector — filtered signal (same length as `x`).
#' @noRd
.lowpass_filter <- function(x, cutoff, sr, order = 8L) {
  x <- as.numeric(x)
  n <- length(x)
  nyq <- sr / 2

  if (cutoff >= nyq) return(x)

  Wn_raw <- cutoff / nyq
  Wn <- min(Wn_raw, 0.99)
  Wn <- max(Wn, 0.01)
  if (!isTRUE(all.equal(Wn_raw, Wn))) {
    warning(sprintf(
      "Normalized cutoff frequency clamped from %.4f to %.4f for numerical stability",
      Wn_raw, Wn
    ), call. = FALSE)
  }

  tryCatch({
    if (!requireNamespace("signal", quietly = TRUE)) {
      stop("signal package not available")
    }
    bf <- signal::butter(order, Wn, type = "low")
    signal::filtfilt(bf, x)
  }, error = function(e) {
    # Fallback: FFT-based brick-wall lowpass
    fft_x <- stats::fft(x)
    freqs <- seq(0, sr - sr / n, length.out = n)
    mask <- freqs <= cutoff | freqs >= (sr - cutoff)
    fft_x[!mask] <- 0
    Re(stats::fft(fft_x, inverse = TRUE)) / n
  })
}

# ---- Signal Pair Extraction -------------------------------------------------

#' Extract a pair of signals from various input types
#'
#' Unifies the way coupling functions receive signal pairs. Accepts:
#' 1. Two numeric vectors (with `sr` required).
#' 2. Two `PhysioExperiment` objects (extracts a single channel from each).
#' 3. A single `MultiPhysioExperiment` with named modalities.
#'
#' When the two signals have different sampling rates (e.g. from different
#' modalities), the function resamples `y` to match `x`'s sampling rate via
#' linear interpolation.
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#' @param y Numeric vector or PhysioExperiment, or NULL when `x` is an MPE.
#' @param sr Numeric — sampling rate (required when x/y are numeric vectors).
#' @param modality_x Character — modality name in MPE for x signal.
#' @param modality_y Character — modality name in MPE for y signal.
#' @param channels_x Integer — which channel to extract from x (default 1).
#' @param channels_y Integer — which channel to extract from y (default 1).
#' @return Named list: `list(x = numeric, y = numeric, sr = numeric)`.
#' @noRd
.extract_signal_pair <- function(x, y = NULL, sr = NULL,
                                 modality_x = NULL, modality_y = NULL,
                                 channels_x = 1L, channels_y = 1L) {

  # Case 1: Two numeric vectors
  if (is.numeric(x) && !is.null(y) && is.numeric(y)) {
    if (is.null(sr)) stop("'sr' is required when x and y are numeric vectors")
    return(list(x = as.numeric(x), y = as.numeric(y), sr = as.numeric(sr)))
  }

  # Case 2: Two PhysioExperiment objects
  if (inherits(x, "PhysioExperiment") && !is.null(y) &&
      inherits(y, "PhysioExperiment")) {
    sr_x <- samplingRate(x)
    sr_y <- samplingRate(y)
    # Extract channel from default (first) assay
    assay_x <- SummarizedExperiment::assay(x, 1L)
    assay_y <- SummarizedExperiment::assay(y, 1L)
    sig_x <- as.numeric(assay_x[, channels_x])
    sig_y <- as.numeric(assay_y[, channels_y])

    # Resample y to match x's sampling rate if different
    if (!is.na(sr_x) && !is.na(sr_y) && abs(sr_x - sr_y) > 1e-6) {
      # Anti-aliasing filter before downsampling
      if (sr_x < sr_y) {
        aa_cutoff <- 0.9 * sr_x / 2
        sig_y <- .lowpass_filter(sig_y, cutoff = aa_cutoff, sr = sr_y)
      }
      sig_y <- .resample_linear(sig_y, sr_from = sr_y, sr_to = sr_x)
      # Trim to same length (take the shorter)
      min_len <- min(length(sig_x), length(sig_y))
      sig_x <- sig_x[seq_len(min_len)]
      sig_y <- sig_y[seq_len(min_len)]
    }

    return(list(x = sig_x, y = sig_y, sr = sr_x))
  }

  # Case 3: MultiPhysioExperiment with modality names
  if (inherits(x, "MultiPhysioExperiment")) {
    if (is.null(modality_x) || is.null(modality_y)) {
      stop("modality_x and modality_y must be specified for MultiPhysioExperiment")
    }

    exps <- experiments(x)
    if (!modality_x %in% names(exps)) {
      stop("modality_x '", modality_x, "' not found in MultiPhysioExperiment")
    }
    if (!modality_y %in% names(exps)) {
      stop("modality_y '", modality_y, "' not found in MultiPhysioExperiment")
    }

    pe_x <- exps[[modality_x]]
    pe_y <- exps[[modality_y]]

    sr_x <- samplingRate(pe_x)
    sr_y <- samplingRate(pe_y)

    assay_x <- SummarizedExperiment::assay(pe_x, 1L)
    assay_y <- SummarizedExperiment::assay(pe_y, 1L)
    sig_x <- as.numeric(assay_x[, channels_x])
    sig_y <- as.numeric(assay_y[, channels_y])

    # The container knows when each stream started. Comparing the head of one
    # array with the head of the other is only correct when they started
    # together; otherwise the two signals are read at different moments and the
    # coupling estimate is of two different time spans.
    off <- .pair_offsets(x, c(modality_x, modality_y))
    t_x <- .pair_times(x, modality_x, sig_x, sr_x, off[[1]])
    t_y <- .pair_times(x, modality_y, sig_y, sr_y, off[[2]])

    # An equal start offset is NOT an equal time correspondence. If either
    # stream carries measured sample times it may have gaps or drift, so the
    # k-th sample of one is not the k-th moment of the other even when both
    # begin together.
    #
    # The historical path is taken only when both streams are on a RECONSTRUCTED
    # REGULAR grid and begin at the same instant. That does not mean the two
    # time vectors are identical -- at 100 Hz and 200 Hz they plainly are not --
    # it means each stream's times are fully determined by its own start and
    # rate, so resampling one onto the other's rate from a shared origin lands
    # on the right moments. This is also the case whose published numbers must
    # not move.
    measured <- .pair_has_measured(x, c(modality_x, modality_y))
    same_start <- is.finite(off[[1]]) && is.finite(off[[2]]) &&
      abs(off[[1]] - off[[2]]) <= 1e-9
    regular_common_start <- !any(measured) && same_start
    span <- NULL

    if (!regular_common_start) {
      lo <- max(t_x[1], t_y[1])
      hi <- min(t_x[length(t_x)], t_y[length(t_y)])
      if (!is.finite(lo) || !is.finite(hi) || hi < lo) {
        stop(sprintf(paste0("'%s' and '%s' do not overlap in time on the shared ",
                            "clock (%s spans [%g, %g] s, %s spans [%g, %g] s). ",
                            "There is no interval on which a coupling estimate ",
                            "would be defined."),
                     modality_x, modality_y,
                     modality_x, t_x[1], t_x[length(t_x)],
                     modality_y, t_y[1], t_y[length(t_y)]), call. = FALSE)
      }
      # Interpolating across a gap would manufacture samples for a stretch the
      # recording does not cover, so a gap inside the window is refused rather
      # than smoothed over.
      .reject_gap_in_window(t_x, sr_x, modality_x, lo, hi)
      .reject_gap_in_window(t_y, sr_y, modality_y, lo, hi)

      grid <- seq(lo, hi, by = 1 / sr_x)
      if (sr_x < sr_y) {
        sig_y <- .lowpass_filter(sig_y, cutoff = 0.9 * sr_x / 2, sr = sr_y)
      }
      sig_x <- stats::approx(t_x, sig_x, xout = grid, method = "linear", rule = 2)$y
      sig_y <- stats::approx(t_y, sig_y, xout = grid, method = "linear", rule = 2)$y
      return(list(x = sig_x, y = sig_y, sr = sr_x,
                  alignment = list(
                    modality_x = modality_x, modality_y = modality_y,
                    offsets = c(off[[1]], off[[2]]),
                    simultaneous_start = same_start,
                    regular_common_start = FALSE,
                    measured_times = measured,
                    overlap = c(lo, hi),
                    output_rate = sr_x,
                    resampled = TRUE,
                    grid = c(first = grid[1], last = grid[length(grid)]),
                    method = paste0("both signals interpolated (linear) onto one ",
                                    "grid over the overlap at the x rate",
                                    if (sr_x < sr_y)
                                      "; y low-pass filtered first" else ""),
                    n_samples = length(grid))))
    }
    simultaneous <- TRUE

    # Resample y to match x if sampling rates differ
    resampled <- FALSE
    if (!is.na(sr_x) && !is.na(sr_y) && abs(sr_x - sr_y) > 1e-6) {
      # Anti-aliasing filter before downsampling
      if (sr_x < sr_y) {
        aa_cutoff <- 0.9 * sr_x / 2
        sig_y <- .lowpass_filter(sig_y, cutoff = aa_cutoff, sr = sr_y)
      }
      sig_y <- .resample_linear(sig_y, sr_from = sr_y, sr_to = sr_x)
      resampled <- TRUE
    }
    min_len <- min(length(sig_x), length(sig_y))
    sig_x <- sig_x[seq_len(min_len)]
    sig_y <- sig_y[seq_len(min_len)]

    return(list(x = sig_x, y = sig_y, sr = sr_x,
                alignment = list(
                  modality_x = modality_x, modality_y = modality_y,
                  offsets = c(off[[1]], off[[2]]),
                  simultaneous_start = simultaneous,
                  regular_common_start = TRUE,
                  measured_times = measured,
                  overlap = span,
                  output_rate = sr_x,
                  resampled = resampled,
                  method = if (resampled)
                    "linear interpolation onto the x rate, low-pass filtered first when downsampling"
                  else "none",
                  n_samples = min_len)))
  }

  stop("Unsupported input types. Provide two numeric vectors, ",
       "two PhysioExperiment objects, or a MultiPhysioExperiment.")
}

# ---- Multitaper Spectral Estimation ------------------------------------------

#' Discrete Prolate Spheroidal Sequences (DPSS / Slepian tapers)
#'
#' Computes DPSS tapers by solving the symmetric tridiagonal eigenvalue
#' problem. Returns the top-k eigenvectors (highest concentration).
#'
#' @param n Integer window length.
#' @param nw Numeric time-bandwidth product (default 4).
#' @param k Integer number of tapers (default \code{2 * nw - 1}).
#' @return Matrix \code{[n x k]} of orthogonal tapers.
#' @noRd
.dpss_tapers <- function(n, nw = 4, k = NULL) {
  if (is.null(k)) k <- as.integer(2 * nw - 1)
  n <- as.integer(n)
  k <- min(k, n)

  # Construct tridiagonal matrix for eigenvalue problem
  diag_main <- ((n - 1 - 2 * seq(0, n - 1))^2 * cos(2 * pi * nw / n)) / 4
  diag_off <- seq_len(n - 1) * (n - seq_len(n - 1)) / 2

  mat <- diag(diag_main)
  for (i in seq_len(n - 1)) {
    mat[i, i + 1] <- diag_off[i]
    mat[i + 1, i] <- diag_off[i]
  }

  ev <- eigen(mat, symmetric = TRUE)
  # Return top-k eigenvectors (highest eigenvalues)
  tapers <- ev$vectors[, seq_len(k), drop = FALSE]

  # Fix sign convention (first element positive)
  for (j in seq_len(k)) {
    if (tapers[1, j] < 0) tapers[, j] <- -tapers[, j]
  }
  tapers
}

#' Multitaper PSD estimate
#'
#' @param x Numeric vector.
#' @param sr Numeric sampling rate.
#' @param nw Numeric time-bandwidth product.
#' @param k Integer number of tapers.
#' @return List with \code{psd} and \code{frequencies}.
#' @noRd
.multitaper_psd <- function(x, sr, nw = 4, k = NULL) {
  n <- length(x)
  tapers <- .dpss_tapers(n, nw, k)
  k_actual <- ncol(tapers)
  freqs <- seq(0, sr / 2, length.out = n %/% 2 + 1)
  n_freqs <- length(freqs)

  psd <- numeric(n_freqs)
  for (j in seq_len(k_actual)) {
    tapered <- x * tapers[, j]
    fft_j <- stats::fft(tapered)
    psd <- psd + Mod(fft_j[seq_len(n_freqs)])^2
  }
  psd <- psd / (k_actual * sr)
  list(psd = psd, frequencies = freqs)
}

#' Multitaper CSD estimate
#'
#' @param x,y Numeric vectors.
#' @param sr Numeric sampling rate.
#' @param nw Numeric time-bandwidth product.
#' @param k Integer number of tapers.
#' @return List with \code{csd} (complex), \code{psd_x}, \code{psd_y},
#'   \code{frequencies}.
#' @noRd
.multitaper_csd <- function(x, y, sr, nw = 4, k = NULL) {
  n <- min(length(x), length(y))
  x <- x[seq_len(n)]
  y <- y[seq_len(n)]
  tapers <- .dpss_tapers(n, nw, k)
  k_actual <- ncol(tapers)
  freqs <- seq(0, sr / 2, length.out = n %/% 2 + 1)
  n_freqs <- length(freqs)

  csd <- complex(n_freqs)
  psd_x <- numeric(n_freqs)
  psd_y <- numeric(n_freqs)

  for (j in seq_len(k_actual)) {
    fx <- stats::fft(x * tapers[, j])
    fy <- stats::fft(y * tapers[, j])
    fx_half <- fx[seq_len(n_freqs)]
    fy_half <- fy[seq_len(n_freqs)]
    csd <- csd + fx_half * Conj(fy_half)
    psd_x <- psd_x + Mod(fx_half)^2
    psd_y <- psd_y + Mod(fy_half)^2
  }

  list(csd = csd / (k_actual * sr),
       psd_x = psd_x / (k_actual * sr),
       psd_y = psd_y / (k_actual * sr),
       frequencies = freqs)
}

# ---- Internal Helpers -------------------------------------------------------

#' Linear interpolation resampling
#'
#' Resamples a signal from one sampling rate to another using linear
#' interpolation (via `stats::approx`).
#'
#' @param x Numeric vector — signal to resample.
#' @param sr_from Numeric — original sampling rate.
#' @param sr_to Numeric — target sampling rate.
#' @return Numeric vector resampled to `sr_to`.
#' @noRd
# Start offsets of two streams on the container's shared clock.
#
# A container saved before the classes were unified has no clock; treating that
# as "both start at zero" would be a guess, so it is reported as NA and the
# caller falls back to the historical head-to-head comparison rather than
# silently claiming an alignment it does not have.
# True sample times of one modality on the shared clock: measured when the
# stream carries measurements, reconstructed from the offset and rate otherwise.
# Does either modality carry measured sample times?
.pair_has_measured <- function(x, modalities) {
  vapply(modalities, function(m)
    isTRUE(tryCatch(PhysioExperiment::hasMeasuredTimes(x, m),
                    error = function(e) FALSE)),
    logical(1), USE.NAMES = TRUE)
}

# Refuse to interpolate across a gap.
#
# THE RULE, AND ITS LIMITS. A step between consecutive measured sample times
# larger than `.GAP_FACTOR` (1.5) times the nominal period 1/rate is treated as a
# gap, and a gap inside the requested window is refused: `approx()` would draw a
# straight line across it and the estimate would rest on samples that were never
# measured.
#
# 1.5 is an ADOPTED CRITERION, not a proof of absence. It sits above the jitter
# of a nominally regular acquisition and below a dropped sample (which gives 2x),
# so it separates those two cases. It does NOT identify every gap under every
# acquisition condition, and two boundaries are worth stating plainly:
#
#   * A stream whose true rate is lower than the rate recorded in its metadata
#     has every step above the threshold, and is refused wholesale rather than
#     analysed. That is the intended direction of failure, but the message will
#     talk about a gap when the real fault is the declared rate.
#   * Drift or irregular sampling that stays within 1.5x is NOT flagged. The
#     window is then interpolated onto the common grid from the real times,
#     which is correct for drift; it is not a claim that the sampling was
#     regular.
#
# For a genuinely irregular stream, resample it explicitly onto a grid you have
# chosen before asking for a coupling estimate.
.GAP_FACTOR <- 1.5

# DETECTION IS NOT DUPLICATED HERE. It is PhysioExperiment::detectGaps(), the
# same function PhysioExperiment::streamCoverage() judges completeness with.
# Sharing the factor was not enough: this guard detected on the whole series and
# then asked which gaps met the window, while streamCoverage() selected the
# window's samples first and differenced those, which loses a gap straddling an
# edge -- its outer bounding sample is outside the window, so the step across it
# is never formed. The two answered differently for the same recording. One
# implementation, one answer.
.reject_gap_in_window <- function(times, rate, modality, lo, hi,
                                  factor = .GAP_FACTOR) {
  g <- PhysioExperiment::detectGaps(times, rate, gap_factor = factor,
                                    window = c(lo, hi))
  if (!nrow(g)) return(invisible(TRUE))
  stop(sprintf(paste0("'%s' has no samples between %g and %g s, inside the ",
                      "requested window [%g, %g]. Interpolating across that gap ",
                      "would invent %.0f samples that were never recorded. ",
                      "Choose a window that avoids the gap, or resample the ",
                      "stream explicitly first."),
               modality, g$start[1L], g$end[1L], lo, hi,
               max(0, (g$end[1L] - g$start[1L]) * rate - 1)), call. = FALSE)
}

.pair_times <- function(x, modality, sig, sr, offset) {
  tt <- tryCatch(PhysioExperiment::streamTimeIndex(x, modality), error = function(e) NULL)
  if (!is.null(tt) && length(tt) == length(sig)) return(as.numeric(tt))
  off <- if (is.finite(offset)) offset else 0
  off + (seq_along(sig) - 1) / sr
}

.pair_offsets <- function(x, modalities) {
  cl <- tryCatch(PhysioExperiment::commonClock(x), error = function(e) NULL)
  off <- cl$offsets
  vapply(modalities, function(m) {
    if (is.null(off) || !m %in% names(off)) NA_real_ else as.numeric(off[[m]])
  }, numeric(1), USE.NAMES = FALSE)
}

.resample_linear <- function(x, sr_from, sr_to) {
  n_from <- length(x)
  duration <- (n_from - 1L) / sr_from
  n_to <- as.integer(round(duration * sr_to)) + 1L

  t_from <- seq(0, duration, length.out = n_from)
  t_to <- seq(0, duration, length.out = n_to)

  stats::approx(t_from, x, xout = t_to, method = "linear",
                rule = 2)$y
}

# What an estimate was computed from. Kept alongside the numeric result so the
# alignment that `.extract_signal_pair()` applied does not stay an internal
# detail: without it a coupling value cannot be related back to the streams,
# channels, time span and conversion that produced it.
.coupling_input_record <- function(pair, sr, channels_x, channels_y, n_samples) {
  al <- pair$alignment
  list(
    modality_x = al$modality_x, modality_y = al$modality_y,
    channels_x = channels_x, channels_y = channels_y,
    offsets = al$offsets,
    simultaneous_start = al$simultaneous_start %||% NA,
    time_span = al$overlap,
    output_rate = sr,
    resampled = al$resampled %||% NA,
    conversion = al$method %||% "none",
    n_samples = n_samples,
    implementation = as.character(utils::packageVersion("PhysioCrossModal"))
  )
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# A stream's measured sample times, or NULL when it is on a regular grid. Read
# from rowData directly so this package does not depend on the foundation's
# internal helper, and kept in step with the column names the foundation treats
# as authoritative.
.MEASURED_TIME_COLUMNS <- c("time_from_t0", "xdf_time", "time_seconds")

.measured_time_column <- function(pe) {
  rd <- tryCatch(SummarizedExperiment::rowData(pe), error = function(e) NULL)
  if (is.null(rd) || nrow(rd) == 0L) return(NULL)
  for (nm in .MEASURED_TIME_COLUMNS) {
    if (nm %in% colnames(rd)) {
      tt <- as.numeric(rd[[nm]])
      if (length(tt) == nrow(rd) && any(is.finite(tt))) return(tt)
    }
  }
  NULL
}
