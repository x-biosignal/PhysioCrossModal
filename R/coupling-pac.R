# Phase-amplitude coupling (PAC): Tort modulation index, Canolty mean vector
# length, Ozkurt normalized MVL, PLV, and the comodulogram. Phase and amplitude
# are extracted with the shared bandpass + Hilbert utilities.

utils::globalVariables(c("phase_freq", "amp_freq", "pac", "display_alpha"))

# Phase of sig_x in phase_band and amplitude envelope of sig_y in amp_band.
.pac_extract <- function(sig_x, sig_y, sr, phase_band, amp_band, order) {
  suppressWarnings({
    phase <- .hilbert_phase(.bandpass_filter(sig_x, phase_band[1], phase_band[2],
                                             sr, order = order))
    amp <- .hilbert_envelope(.bandpass_filter(sig_y, amp_band[1], amp_band[2],
                                              sr, order = order))
  })
  list(phase = phase, amp = amp)
}

# Tort Kullback-Leibler modulation index over phase bins (in 0..1).
.pac_tort <- function(phase, amp, n_bins = 18L) {
  br <- seq(-pi, pi, length.out = n_bins + 1L)
  bi <- cut(phase, br, include.lowest = TRUE, labels = FALSE)
  ma <- tapply(amp, factor(bi, levels = seq_len(n_bins)), mean)
  ma[is.na(ma)] <- 0
  total <- sum(ma)
  if (total <= 0) return(list(mi = 0, dist = rep(1 / n_bins, n_bins)))
  p <- ma / total
  positive <- p > 0
  mi <- (log(n_bins) + sum(p[positive] * log(p[positive]))) / log(n_bins)
  list(mi = as.numeric(mi), dist = as.numeric(p))
}

# Canolty mean vector length (raw, amplitude-scale dependent).
.pac_canolty <- function(phase, amp) Mod(mean(amp * exp(1i * phase)))

# Ozkurt normalized mean vector length (in 0..1).
.pac_ozkurt <- function(phase, amp) {
  n <- length(amp)
  denom <- sqrt(n) * sqrt(sum(amp^2))
  if (denom <= 0) return(0)
  Mod(sum(amp * exp(1i * phase))) / denom
}

# Phase-locking value between the low-frequency phase and the phase of the
# amplitude envelope (band-passed at the phase band).
.pac_plv <- function(phase, amp, sr, phase_band, order) {
  suppressWarnings(
    amp_phase <- .hilbert_phase(.bandpass_filter(amp, phase_band[1],
                                                 phase_band[2], sr, order = order)))
  Mod(mean(exp(1i * (phase - amp_phase))))
}

.pac_stat <- function(phase, amp, method, n_bins, sr, phase_band, order) {
  switch(method,
    tort = .pac_tort(phase, amp, n_bins)$mi,
    canolty = .pac_canolty(phase, amp),
    ozkurt = .pac_ozkurt(phase, amp),
    plv = .pac_plv(phase, amp, sr, phase_band, order))
}

.pac_positive_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop("'", name, "' must be a finite integer >= ", minimum, ".",
         call. = FALSE)
  }
  as.integer(x)
}

.pac_validate_frequency_vector <- function(x, name) {
  if (!is.numeric(x) || length(x) < 1L || any(!is.finite(x)) ||
      anyDuplicated(x) || (length(x) > 1L && any(diff(x) <= 0))) {
    stop("'", name, "' must be a non-empty, finite, strictly increasing, ",
         "unique numeric vector.", call. = FALSE)
  }
  as.numeric(x)
}

.pac_realized_bands <- function(freqs, bandwidth, nyquist, name) {
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
      !is.finite(bandwidth) || bandwidth <= 0) {
    stop("'", name, "_bw' must be one positive finite number.", call. = FALSE)
  }
  requested_low <- freqs - bandwidth
  requested_high <- freqs + bandwidth
  if (any(requested_low <= 0) || any(requested_high >= nyquist)) {
    stop("Every realized ", name,
         " band must lie strictly within (0, Nyquist).", call. = FALSE)
  }
  bands <- cbind(
    center = freqs,
    requested_low = requested_low,
    requested_high = requested_high,
    low = pmax(requested_low, 0.01 * nyquist),
    high = pmin(requested_high, 0.99 * nyquist)
  )
  rownames(bands) <- as.character(freqs)
  if (any(bands[, "low"] >= bands[, "high"])) {
    stop("A realized ", name, " band collapsed after numerical boundary ",
         "handling.", call. = FALSE)
  }
  if (anyDuplicated(data.frame(low = bands[, "low"],
                               high = bands[, "high"]))) {
    stop("Realized ", name, " bands must be unique.", call. = FALSE)
  }

  bands
}

.pac_filter_phases <- function(signal, bands, sr, order) {
  suppressWarnings(lapply(seq_len(nrow(bands)), function(i) {
    .hilbert_phase(.bandpass_filter(signal, bands[i, "low"],
                                    bands[i, "high"], sr, order,
                                    clamp = FALSE))
  }))
}

.pac_filter_amplitudes <- function(signal, bands, sr, order) {
  suppressWarnings(lapply(seq_len(nrow(bands)), function(i) {
    .hilbert_envelope(.bandpass_filter(signal, bands[i, "low"],
                                       bands[i, "high"], sr, order,
                                       clamp = FALSE))
  }))
}

.pac_grid_from_components <- function(phases, amplitudes, method, n_bins, sr,
                                      phase_bands, order) {
  phase_names <- rownames(phase_bands)
  amp_names <- names(amplitudes)
  mat <- matrix(NA_real_, length(phases), length(amplitudes),
                dimnames = list(phase_names, amp_names))
  for (i in seq_along(phases)) {
    for (j in seq_along(amplitudes)) {
      value <- .pac_stat(
        phases[[i]], amplitudes[[j]], method, n_bins, sr,
        phase_bands[i, c("low", "high")], order
      )
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
        stop("PAC computation produced a non-finite scalar.", call. = FALSE)
      }
      if (method %in% c("tort", "ozkurt")) {
        tolerance <- 64 * .Machine$double.eps
        if (value < -tolerance || value > 1 + tolerance) {
          stop("Normalized PAC fell outside its theoretical [0, 1] range.",
               call. = FALSE)
        }
        value <- min(1, max(0, value))
      } else if (value < 0) {
        stop("PAC computation produced a negative value.", call. = FALSE)
      }
      mat[i, j] <- value
    }
  }
  mat
}

.pac_observed_grid <- function(sig_x, sig_y, sr, phase_bands, amp_bands,
                               method, n_bins, order) {
  phases <- .pac_filter_phases(sig_x, phase_bands, sr, order)
  amplitudes <- .pac_filter_amplitudes(sig_y, amp_bands, sr, order)
  names(amplitudes) <- rownames(amp_bands)
  list(
    matrix = .pac_grid_from_components(phases, amplitudes, method, n_bins, sr,
                                       phase_bands, order),
    phases = phases
  )
}

.pac_peak <- function(mat, phase_freqs, amp_freqs) {
  pk <- which(mat == max(mat, na.rm = TRUE), arr.ind = TRUE)[1L, ]
  list(phase_freq = phase_freqs[pk[1L]], amp_freq = amp_freqs[pk[2L]],
       value = mat[pk[1L], pk[2L]])
}

.pac_restore_rng <- function(had_seed, old_seed) {
  if (had_seed) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

#' Phase-amplitude coupling
#'
#' Quantifies how the amplitude of a fast oscillation (in \code{amp_band}) is
#' modulated by the phase of a slow oscillation (in \code{phase_band}). The
#' phase is taken from \code{x} and the amplitude envelope from \code{y} (or
#' from \code{x} itself when \code{y} is \code{NULL}), each via a bandpass filter
#' and the Hilbert transform. Supported measures are the Tort modulation index
#' (a normalized Kullback-Leibler distance of the phase-amplitude distribution
#' from uniform, in 0..1), the Canolty mean vector length (raw), the Ozkurt
#' normalized mean vector length (in 0..1), and a phase-locking value.
#'
#' @param x A numeric signal or PhysioExperiment (phase source), or a
#'   MultiPhysioExperiment.
#' @param y A numeric signal or PhysioExperiment (amplitude source), or
#'   \code{NULL} to use \code{x} for both.
#' @param sr Sampling rate in Hz (required for numeric input).
#' @param phase_band Numeric length-2 phase frequency band in Hz
#'   (default: \code{c(4, 8)}).
#' @param amp_band Numeric length-2 amplitude frequency band in Hz
#'   (default: \code{c(30, 80)}).
#' @param method \code{"tort"}, \code{"canolty"}, \code{"ozkurt"}, or
#'   \code{"plv"}.
#' @param n_bins Number of phase bins for the Tort modulation index
#'   (default: 18).
#' @param order Bandpass filter order (default: 4).
#' @param modality_x,modality_y,channels_x,channels_y Passed to the shared
#'   signal-pair extraction.
#' @return A list with \code{pac} (the coupling value), \code{method},
#'   \code{phase_band}, \code{amp_band}, and, for the Tort method, the
#'   per-bin normalized amplitude \code{distribution} and \code{bin_centers}.
#' @references
#' Tort, A. B. L., et al. (2010). Measuring phase-amplitude coupling between
#' neuronal oscillations of different frequencies. Journal of Neurophysiology,
#' 104(2), 1195-1210.
#'
#' Canolty, R. T., et al. (2006). High gamma power is phase-locked to theta
#' oscillations in human neocortex. Science, 313(5793), 1626-1628.
#' @seealso [comodulogram()], [plotComodulogram()], [surrogateTest()]
#' @export
#' @examples
#' set.seed(1)
#' t <- (0:4999) / 500
#' ph <- 2 * pi * 6 * t
#' sig <- sin(ph) + (0.1 + (1 + cos(ph - pi)) / 2) * sin(2 * pi * 50 * t)
#' phaseAmplitudeCoupling(sig, sr = 500, phase_band = c(4, 8),
#'                        amp_band = c(35, 65))$pac
phaseAmplitudeCoupling <- function(x, y = NULL, sr = NULL,
                                   phase_band = c(4, 8), amp_band = c(30, 80),
                                   method = c("tort", "canolty", "ozkurt", "plv"),
                                   n_bins = 18L, order = 4L,
                                   modality_x = NULL, modality_y = NULL,
                                   channels_x = 1L, channels_y = 1L) {
  method <- match.arg(method)
  stopifnot(is.numeric(phase_band), length(phase_band) == 2,
            phase_band[1] < phase_band[2])
  stopifnot(is.numeric(amp_band), length(amp_band) == 2, amp_band[1] < amp_band[2])
  if (is.null(y) && is.numeric(x)) y <- x    # single-signal PAC: phase and amp from x
  pair <- .extract_signal_pair(x, y, sr, modality_x = modality_x,
                               modality_y = modality_y,
                               channels_x = channels_x, channels_y = channels_y)
  pa <- .pac_extract(pair$x, pair$y, pair$sr, phase_band, amp_band, order)
  out <- list(pac = NA_real_, method = method, phase_band = phase_band,
              amp_band = amp_band, n_bins = as.integer(n_bins))
  if (method == "tort") {
    tt <- .pac_tort(pa$phase, pa$amp, n_bins)
    out$pac <- tt$mi
    out$distribution <- tt$dist
    br <- seq(-pi, pi, length.out = n_bins + 1L)
    out$bin_centers <- (br[-1] + br[-(n_bins + 1L)]) / 2
  } else {
    out$pac <- .pac_stat(pa$phase, pa$amp, method, n_bins, pair$sr, phase_band, order)
  }
  out
}

#' Comodulogram (phase-frequency by amplitude-frequency PAC)
#'
#' Computes phase-amplitude coupling over a grid of phase frequencies and
#' amplitude frequencies, yielding a comodulogram matrix whose peak locates the
#' dominant coupling.
#'
#' @param x,y,sr,method,n_bins,order,modality_x,modality_y,channels_x,channels_y
#'   As in [phaseAmplitudeCoupling()].
#' @param phase_freqs Numeric vector of phase centre frequencies in Hz.
#' @param amp_freqs Numeric vector of amplitude centre frequencies in Hz.
#' @param phase_bw Half-bandwidth of each phase band in Hz (default: 2).
#' @param amp_bw Half-bandwidth of each amplitude band in Hz (default: 10).
#' @return A list with the \code{matrix} (rows = phase frequencies, columns =
#'   amplitude frequencies), \code{phase_freqs}, \code{amp_freqs},
#'   \code{method}, and \code{peak} (the (phase, amplitude) frequency of the
#'   maximum and its value).
#' @seealso [phaseAmplitudeCoupling()], [plotComodulogram()]
#' @export
#' @examples
#' \dontrun{
#' set.seed(1); t <- (0:4999) / 500; ph <- 2 * pi * 6 * t
#' sig <- sin(ph) + (0.1 + (1 + cos(ph - pi)) / 2) * sin(2 * pi * 50 * t)
#' cm <- comodulogram(sig, sr = 500, phase_freqs = c(4, 6, 8),
#'                    amp_freqs = c(30, 50, 70))
#' cm$peak
#' }
comodulogram <- function(x, y = NULL, sr = NULL, phase_freqs, amp_freqs,
                         method = c("tort", "canolty", "ozkurt", "plv"),
                         phase_bw = 2, amp_bw = 10, n_bins = 18L, order = 4L,
                         modality_x = NULL, modality_y = NULL,
                         channels_x = 1L, channels_y = 1L) {
  method <- match.arg(method)
  stopifnot(is.numeric(phase_freqs), is.numeric(amp_freqs))
  if (is.null(y) && is.numeric(x)) y <- x
  pair <- .extract_signal_pair(x, y, sr, modality_x = modality_x,
                               modality_y = modality_y,
                               channels_x = channels_x, channels_y = channels_y)
  nyq <- pair$sr / 2
  phase_bands <- cbind(
    center = phase_freqs,
    low = pmax(0.5, 0.01 * nyq, phase_freqs - phase_bw),
    high = pmin(nyq - 0.5, 0.99 * nyq, phase_freqs + phase_bw)
  )
  amp_bands <- cbind(
    center = amp_freqs,
    low = pmax(0.5, 0.01 * nyq, amp_freqs - amp_bw),
    high = pmin(nyq - 0.5, 0.99 * nyq, amp_freqs + amp_bw)
  )
  rownames(phase_bands) <- as.character(phase_freqs)
  rownames(amp_bands) <- as.character(amp_freqs)
  observed <- .pac_observed_grid(pair$x, pair$y, pair$sr, phase_bands,
                                 amp_bands, method, n_bins, order)
  mat <- observed$matrix
  list(matrix = mat, phase_freqs = phase_freqs, amp_freqs = amp_freqs,
       method = method, peak = .pac_peak(mat, phase_freqs, amp_freqs))
}

#' Modulation-index comodulogram with surrogate inference
#'
#' Computes Tort, Canolty, or Ozkurt phase-amplitude coupling over a declared
#' frequency grid. Each surrogate amplitude-source signal produces one complete
#' matrix, so cell-wise p-values share a coherent null replicate. P-values use
#' the conservative Phipson-Smyth correction and are adjusted over the complete
#' grid.
#'
#' PAC is an association and can be produced by waveform shape, filter leakage,
#' edge transients, non-stationarity, or common inputs. A high or significant
#' value does not establish cross-frequency communication or a biological
#' mechanism. The returned threshold is the unadjusted cell-wise null quantile;
#' the \code{significant} mask is instead defined from adjusted p-values.
#' Benjamini-Hochberg controls the false discovery rate, not strong family-wise
#' error. Requested and exact realized pass bands are both retained; a
#' structured diagnostic records any adjustment at the filter's numerical
#' boundary.
#'
#' @param x A numeric signal or PhysioExperiment providing the phase source, or
#'   a MultiPhysioExperiment.
#' @param phase_freqs Strictly increasing phase-band centre frequencies in Hz.
#' @param amp_freqs Strictly increasing amplitude-band centre frequencies in Hz.
#' @param sr Sampling rate in Hz, required for numeric input.
#' @param y Optional amplitude-source signal or PhysioExperiment. When
#'   \code{NULL}, a non-MultiPhysioExperiment \code{x} supplies both sources.
#' @param method \code{"tort"}, \code{"canolty"}, or \code{"ozkurt"}.
#' @param phase_bw,amp_bw Half-bandwidths in Hz.
#' @param n_bins Number of phase bins for Tort MI.
#' @param order Bandpass filter order.
#' @param n_surrogates Number of full-grid surrogates. Use zero to disable
#'   inference explicitly.
#' @param surrogate_type \code{"timeshift"}, \code{"phase"}, \code{"iaaft"},
#'   or \code{"aaft"}.
#' @param alpha Significance level.
#' @param p_adjust_method A method from \code{stats::p.adjust.methods}, applied
#'   once across the complete frequency grid.
#' @param seed Optional non-negative integer seed.
#' @param cores Positive integer number of worker processes.
#' @param modality_x,modality_y,channels_x,channels_y Passed to the shared
#'   signal-pair extractor.
#' @return A \code{pac_comodulogram} list containing the observed matrix,
#'   grid and realized bands, peak, cell-wise p-values, adjusted p-values,
#'   significance mask, null thresholds and summaries, surrogate settings, and
#'   structured diagnostic strings. Inferential matrices are \code{NULL} when
#'   \code{n_surrogates = 0}.
#' @references
#' Tort, A. B. L., et al. (2010). Journal of Neurophysiology, 104, 1195-1210.
#'
#' Canolty, R. T., et al. (2006). Science, 313, 1626-1628.
#'
#' Ozkurt, T. E., & Schnitzler, A. (2011). Journal of Neuroscience Methods,
#' 201, 438-443.
#'
#' Phipson, B., & Smyth, G. K. (2010). Statistical Applications in Genetics
#' and Molecular Biology, 9, Article 39.
#' @seealso [comodulogram()], [phaseAmplitudeCoupling()],
#'   [plotComodulogram()]
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' t <- (0:4999) / 500
#' ph <- 2 * pi * 6 * t
#' sig <- sin(ph) + (0.1 + (1 + cos(ph - pi)) / 2) * sin(2 * pi * 50 * t)
#' mi <- modulationIndex(sig, c(4, 6, 8), c(40, 50, 60), sr = 500,
#'                       n_surrogates = 99, seed = 1)
#' plotComodulogram(mi)
#' }
modulationIndex <- function(
    x, phase_freqs, amp_freqs, sr = NULL, y = NULL,
    method = c("tort", "canolty", "ozkurt"),
    phase_bw = 2, amp_bw = 10, n_bins = 18L, order = 4L,
    n_surrogates = 199L,
    surrogate_type = c("timeshift", "phase", "iaaft", "aaft"),
    alpha = 0.05, p_adjust_method = "BH", seed = NULL, cores = 1L,
    modality_x = NULL, modality_y = NULL,
    channels_x = 1L, channels_y = 1L) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit(.pac_restore_rng(had_seed, old_seed), add = TRUE)

  method <- match.arg(method)
  surrogate_type <- match.arg(surrogate_type)
  phase_freqs <- .pac_validate_frequency_vector(phase_freqs, "phase_freqs")
  amp_freqs <- .pac_validate_frequency_vector(amp_freqs, "amp_freqs")
  n_bins <- .pac_positive_integer(n_bins, "n_bins", minimum = 2L)
  order <- .pac_positive_integer(order, "order")
  n_surrogates <- .pac_positive_integer(n_surrogates, "n_surrogates",
                                        minimum = 0L)
  cores <- .pac_positive_integer(cores, "cores")
  channels_x <- .pac_positive_integer(channels_x, "channels_x")
  channels_y <- .pac_positive_integer(channels_y, "channels_y")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) ||
      !p_adjust_method %in% stats::p.adjust.methods) {
    stop("'p_adjust_method' must be one of stats::p.adjust.methods.",
         call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .pac_positive_integer(seed, "seed", minimum = 0L)
  }
  minimum_surrogates <- ceiling(1 / alpha)
  if (n_surrogates > 0L && n_surrogates < minimum_surrogates) {
    stop("'n_surrogates' must be zero or at least ", minimum_surrogates,
         " for the requested alpha.", call. = FALSE)
  }

  if (is.numeric(x) && !is.null(dim(x))) {
    stop("Numeric 'x' must be a vector, not an array or matrix.", call. = FALSE)
  }
  if (!is.null(y) && is.numeric(y) && !is.null(dim(y))) {
    stop("Numeric 'y' must be a vector, not an array or matrix.", call. = FALSE)
  }
  if (is.null(y) && !inherits(x, "MultiPhysioExperiment")) y <- x
  pair <- .extract_signal_pair(
    x, y, sr, modality_x = modality_x, modality_y = modality_y,
    channels_x = channels_x, channels_y = channels_y
  )
  if (!is.numeric(pair$sr) || length(pair$sr) != 1L ||
      !is.finite(pair$sr) || pair$sr <= 0) {
    stop("The realized sampling rate must be one positive finite number.",
         call. = FALSE)
  }
  if (length(pair$x) != length(pair$y)) {
    stop("Phase and amplitude source signals must have equal lengths.",
         call. = FALSE)
  }
  if (length(pair$x) < 2L || any(!is.finite(pair$x)) ||
      any(!is.finite(pair$y))) {
    stop("Extracted signals must be finite vectors with at least two samples.",
         call. = FALSE)
  }
  if (diff(range(pair$x)) == 0 || diff(range(pair$y)) == 0) {
    stop("Extracted signals must be non-constant.", call. = FALSE)
  }
  min_length <- 3L * (2L * order + 1L) + 1L
  if (length(pair$x) < min_length) {
    stop("Extracted signals are too short for the requested filter order; ",
         "need at least ", min_length, " samples.", call. = FALSE)
  }

  nyquist <- pair$sr / 2
  phase_bands <- .pac_realized_bands(phase_freqs, phase_bw, nyquist, "phase")
  amp_bands <- .pac_realized_bands(amp_freqs, amp_bw, nyquist, "amplitude")
  observed <- .pac_observed_grid(
    pair$x, pair$y, pair$sr, phase_bands, amp_bands, method, n_bins, order
  )
  observed_matrix <- observed$matrix
  warnings <- c(
    "interpretation:pac_is_association_not_mechanism",
    "diagnostic:waveform_filter_edge_and_nonstationarity_checks_required"
  )
  if (any(phase_bands[, "requested_low"] != phase_bands[, "low"]) ||
      any(phase_bands[, "requested_high"] != phase_bands[, "high"]) ||
      any(amp_bands[, "requested_low"] != amp_bands[, "low"]) ||
      any(amp_bands[, "requested_high"] != amp_bands[, "high"])) {
    warnings <- c(warnings,
                  "filter:requested_band_adjusted_at_numerical_boundary")
  }
  overlaps <- outer(
    phase_bands[, "low"], amp_bands[, "high"], "<"
  ) & outer(phase_bands[, "high"], amp_bands[, "low"], ">")
  if (any(overlaps)) {
    warnings <- c(warnings,
                  "diagnostic:phase_and_amplitude_pass_bands_overlap")
  }

  out <- list(
    matrix = observed_matrix,
    phase_freqs = phase_freqs,
    amp_freqs = amp_freqs,
    method = method,
    phase_bw = as.numeric(phase_bw),
    amp_bw = as.numeric(amp_bw),
    peak = .pac_peak(observed_matrix, phase_freqs, amp_freqs),
    p_value = NULL,
    p_adjusted = NULL,
    significant = NULL,
    threshold = NULL,
    surrogate_summary = NULL,
    surrogate_settings = list(
      type = surrogate_type, count = n_surrogates, alpha = alpha,
      adjustment = p_adjust_method, seed = seed, cores = cores,
      quantile_type = 8L
    ),
    realized_bands = list(phase = phase_bands, amplitude = amp_bands),
    warnings = warnings
  )

  if (n_surrogates == 0L) {
    out$warnings <- c(out$warnings,
                      "inference_disabled:n_surrogates_is_zero")
    class(out) <- c("pac_comodulogram", "list")
    return(out)
  }

  if (!is.null(seed)) set.seed(seed)
  generator <- switch(
    surrogate_type,
    timeshift = .timeshift_surrogate,
    phase = .phase_randomize_surrogate,
    iaaft = .iaaft_surrogate,
    aaft = .aaft_surrogate
  )

  # Generate all stochastic inputs in a fixed order before any parallel work.
  surrogate_signals <- lapply(seq_len(n_surrogates), function(i) {
    generator(pair$y)
  })
  surrogate_grids <- .parallel_apply(
    surrogate_signals,
    function(y_surrogate) {
      amplitudes <- .pac_filter_amplitudes(
        y_surrogate, amp_bands, pair$sr, order
      )
      names(amplitudes) <- rownames(amp_bands)
      .pac_grid_from_components(
        observed$phases, amplitudes, method, n_bins, pair$sr,
        phase_bands, order
      )
    },
    cores = cores
  )
  n_cells <- length(observed_matrix)
  null_flat <- matrix(
    vapply(surrogate_grids, as.vector, numeric(n_cells)),
    nrow = n_cells, ncol = n_surrogates
  )
  dimnames_grid <- dimnames(observed_matrix)
  as_grid <- function(values) {
    matrix(values, nrow = nrow(observed_matrix),
           ncol = ncol(observed_matrix), dimnames = dimnames_grid)
  }

  observed_flat <- as.vector(observed_matrix)
  exceedances <- rowSums(null_flat >= observed_flat)
  p_value <- as_grid((exceedances + 1) / (n_surrogates + 1))
  p_adjusted <- as_grid(stats::p.adjust(as.vector(p_value),
                                        method = p_adjust_method))
  null_mean <- as_grid(rowMeans(null_flat))
  null_sd <- as_grid(apply(null_flat, 1L, stats::sd))
  threshold <- as_grid(apply(
    null_flat, 1L, stats::quantile, probs = 1 - alpha,
    names = FALSE, type = 8L
  ))

  out$p_value <- p_value
  out$p_adjusted <- p_adjusted
  out$significant <- p_adjusted <= alpha
  out$threshold <- threshold
  out$surrogate_summary <- list(
    mean = null_mean, sd = null_sd, quantile = threshold
  )
  if (identical(p_adjust_method, "BH")) {
    out$warnings <- c(out$warnings,
                      "multiplicity:BH_controls_FDR_not_familywise_error")
  }
  if (isTRUE(getOption("PhysioCrossModal.pac.keep_surrogates", FALSE))) {
    out$surrogate_array <- array(
      null_flat,
      dim = c(nrow(observed_matrix), ncol(observed_matrix), n_surrogates),
      dimnames = c(dimnames_grid, list(as.character(seq_len(n_surrogates))))
    )
  }
  class(out) <- c("pac_comodulogram", "list")
  out
}

#' Plot a comodulogram
#'
#' Renders a comodulogram (from [comodulogram()]) as a heatmap.
#'
#' @param result A comodulogram result list.
#' @param title Plot title.
#' @param mask Optional logical significance mask or numeric p-value matrix.
#'   When \code{NULL}, uses \code{result$significant} if present.
#' @param nonsignificant_alpha Opacity for non-significant cells.
#' @return A \code{ggplot} object.
#' @seealso [comodulogram()], [modulationIndex()]
#' @export
#' @examples
#' \dontrun{
#' plotComodulogram(comodulogram(sig, sr = 500,
#'   phase_freqs = c(4, 6, 8), amp_freqs = c(30, 50, 70)))
#' }
plotComodulogram <- function(result, title = "Comodulogram", mask = NULL,
                             nonsignificant_alpha = 0.2) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotComodulogram().", call. = FALSE)
  }
  if (!is.list(result) || !is.matrix(result$matrix) ||
      !is.numeric(result$matrix) || any(!is.finite(result$matrix))) {
    stop("'result$matrix' must be a finite numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(result$phase_freqs) || !is.numeric(result$amp_freqs) ||
      length(result$phase_freqs) != nrow(result$matrix) ||
      length(result$amp_freqs) != ncol(result$matrix)) {
    stop("Frequency vectors must match 'result$matrix'.", call. = FALSE)
  }
  if (!is.numeric(nonsignificant_alpha) ||
      length(nonsignificant_alpha) != 1L ||
      !is.finite(nonsignificant_alpha) ||
      nonsignificant_alpha < 0 || nonsignificant_alpha > 1) {
    stop("'nonsignificant_alpha' must be between zero and one.",
         call. = FALSE)
  }
  if (is.null(mask) && !is.null(result$significant)) {
    mask <- result$significant
  }
  if (is.null(mask)) {
    mask <- matrix(TRUE, nrow(result$matrix), ncol(result$matrix),
                   dimnames = dimnames(result$matrix))
  } else {
    if (!is.matrix(mask) || !identical(dim(mask), dim(result$matrix)) ||
        !identical(dimnames(mask), dimnames(result$matrix))) {
      stop("'mask' dimensions and dimnames must match 'result$matrix'.",
           call. = FALSE)
    }
    if (is.logical(mask)) {
      if (anyNA(mask)) stop("Logical 'mask' cannot contain NA.", call. = FALSE)
    } else if (is.numeric(mask)) {
      alpha <- result$surrogate_settings$alpha
      if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
          alpha <= 0 || alpha >= 1) {
        stop("A numeric 'mask' requires a valid result significance alpha.",
             call. = FALSE)
      }
      if (any(!is.finite(mask)) || any(mask < 0 | mask > 1)) {
        stop("A numeric 'mask' must contain finite p-values in [0, 1].",
             call. = FALSE)
      }
      mask <- mask <= alpha
    } else {
      stop("'mask' must be a logical or numeric matrix.", call. = FALSE)
    }
  }

  df <- expand.grid(phase_freq = result$phase_freqs, amp_freq = result$amp_freqs)
  df$pac <- as.vector(result$matrix)
  df$display_alpha <- ifelse(as.vector(mask), 1, nonsignificant_alpha)
  ggplot2::ggplot(df, ggplot2::aes(x = phase_freq, y = amp_freq, fill = pac)) +
    ggplot2::geom_tile(ggplot2::aes(alpha = display_alpha)) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_fill_viridis_c(option = "inferno", name = result$method) +
    ggplot2::labs(x = "Phase frequency (Hz)", y = "Amplitude frequency (Hz)",
                  title = title) +
    ggplot2::theme_minimal()
}
