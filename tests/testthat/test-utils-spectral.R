# Tests for internal spectral utility functions (utils-spectral.R)

# ---- Window Functions -------------------------------------------------------

test_that(".window_hanning produces known values for n=4", {
  w <- PhysioCrossModal:::.window_hanning(4)
  expect_equal(w, c(0, 0.75, 0.75, 0), tolerance = 1e-12)
})

test_that(".window_hanning edge cases", {
  expect_equal(PhysioCrossModal:::.window_hanning(1), 1)
  expect_length(PhysioCrossModal:::.window_hanning(0), 0)
  w5 <- PhysioCrossModal:::.window_hanning(5)
  expect_length(w5, 5)
  expect_equal(w5[1], 0)
  expect_equal(w5[5], 0)
  # Symmetric
  expect_equal(w5[2], w5[4], tolerance = 1e-12)
  expect_equal(w5[3], 1)  # peak at center for odd n
})

test_that(".window_hamming produces correct values for n=4", {
  w <- PhysioCrossModal:::.window_hamming(4)
  expected <- 0.54 - 0.46 * cos(2 * pi * (0:3) / 3)
  expect_equal(w, expected, tolerance = 1e-12)
})

test_that(".window_hamming edge cases", {
  expect_equal(PhysioCrossModal:::.window_hamming(1), 1)
  expect_length(PhysioCrossModal:::.window_hamming(0), 0)
})

# ---- Welch PSD --------------------------------------------------------------

test_that(".welch_psd detects peak at 10Hz for sinusoidal input", {
  sr <- 500
  n_sec <- 4
  n <- sr * n_sec
  t <- seq(0, n_sec - 1 / sr, length.out = n)
  x <- sin(2 * pi * 10 * t)

  result <- PhysioCrossModal:::.welch_psd(x, nperseg = 256L, sr = sr)

  expect_type(result, "list")
  expect_true("psd" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_length(result$psd, length(result$frequencies))

  # Peak frequency should be at 10Hz
  peak_idx <- which.max(result$psd)
  peak_freq <- result$frequencies[peak_idx]
  expect_true(abs(peak_freq - 10) < 2,
    info = paste("Peak at", peak_freq, "Hz, expected ~10 Hz"))
})

test_that(".welch_psd returns positive values", {
  sr <- 100
  x <- rnorm(1000)
  result <- PhysioCrossModal:::.welch_psd(x, nperseg = 64L, sr = sr)
  expect_true(all(result$psd >= 0))
})

test_that(".welch_psd handles short signals by adapting segment length", {
  x <- rnorm(50)
  result <- PhysioCrossModal:::.welch_psd(x, nperseg = 256L, sr = 100)
  expect_type(result, "list")
  expect_true(length(result$psd) > 0)
})

test_that(".welch_psd with hamming window works", {
  x <- sin(2 * pi * 10 * seq(0, 2, length.out = 1000))
  result <- PhysioCrossModal:::.welch_psd(x, nperseg = 128L, sr = 500,
                                           window = "hamming")
  expect_type(result, "list")
  peak_idx <- which.max(result$psd)
  peak_freq <- result$frequencies[peak_idx]
  expect_true(abs(peak_freq - 10) < 4)
})

test_that(".welch_psd validates inputs", {
  expect_error(PhysioCrossModal:::.welch_psd(rnorm(100), nperseg = 1L),
               "nperseg")
  expect_error(
    PhysioCrossModal:::.welch_psd(rnorm(100), nperseg = 64L, noverlap = -1L),
    "noverlap"
  )
})

# ---- Welch CSD --------------------------------------------------------------

test_that(".welch_csd of identical signals is real-valued and matches PSD", {
  sr <- 500
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.1 * rnorm(n)

  psd_result <- PhysioCrossModal:::.welch_psd(x, nperseg = 256L, sr = sr)
  csd_result <- PhysioCrossModal:::.welch_csd(x, x, nperseg = 256L, sr = sr)

  expect_length(csd_result$csd, length(psd_result$psd))
  expect_equal(csd_result$frequencies, psd_result$frequencies)

  # CSD of identical signals should be real-valued (imaginary part ~0)
  imag_ratio <- max(abs(Im(csd_result$csd))) / max(abs(Re(csd_result$csd)))
  expect_true(imag_ratio < 1e-10,
    info = paste("Imaginary/Real ratio:", imag_ratio))

  # Real part of CSD should match PSD (before one-sided doubling)
  # The PSD applies one-sided doubling but CSD does not, so we compare shapes
  # Both should have the peak at the same frequency
  csd_peak <- which.max(Mod(csd_result$csd))
  psd_peak <- which.max(psd_result$psd)
  expect_equal(csd_result$frequencies[csd_peak],
               psd_result$frequencies[psd_peak])
})

test_that(".welch_csd rejects unequal-length signals", {
  expect_error(
    PhysioCrossModal:::.welch_csd(rnorm(100), rnorm(200)),
    "same length"
  )
})

test_that(".welch_csd returns complex vector", {
  x <- rnorm(500)
  y <- rnorm(500)
  result <- PhysioCrossModal:::.welch_csd(x, y, nperseg = 64L, sr = 100)
  expect_true(is.complex(result$csd))
})

# ---- Hilbert Transform ------------------------------------------------------

test_that(".hilbert_analytic of cosine has envelope ~1.0", {
  sr <- 1000
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- cos(2 * pi * 10 * t)

  analytic <- PhysioCrossModal:::.hilbert_analytic(x)

  expect_length(analytic, n)
  expect_true(is.complex(analytic))

  # Real part should equal the original signal
  expect_equal(Re(analytic), x, tolerance = 1e-10)

  # Envelope should be ~1.0 (away from edges to avoid boundary effects)
  envelope <- Mod(analytic)
  # Check central 80% of the signal
  inner <- seq(floor(n * 0.1), floor(n * 0.9))
  expect_true(all(abs(envelope[inner] - 1.0) < 0.05),
    info = paste("Max envelope deviation from 1.0:",
                 max(abs(envelope[inner] - 1.0))))
})

test_that(".hilbert_phase of cosine matches expected phase progression", {
  sr <- 1000
  n <- 1000
  t <- seq(0, (n - 1) / sr, length.out = n)
  freq <- 5
  x <- cos(2 * pi * freq * t)

  phase <- PhysioCrossModal:::.hilbert_phase(x)

  expect_length(phase, n)
  expect_true(is.numeric(phase))

  # Phase should be within [-pi, pi]
  expect_true(all(phase >= -pi - 1e-10 & phase <= pi + 1e-10))

  # At t=0, cos(0) = 1, so the analytic signal is 1 + 0i, phase should be ~0
  expect_true(abs(phase[1]) < 0.1, info = paste("Phase at t=0:", phase[1]))
})

test_that(".hilbert_envelope of cosine is approximately constant", {
  sr <- 500
  n <- 1000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- cos(2 * pi * 20 * t)

  envelope <- PhysioCrossModal:::.hilbert_envelope(x)
  expect_length(envelope, n)

  # Central 80% should have envelope ~ 1.0
  inner <- seq(floor(n * 0.1), floor(n * 0.9))
  expect_true(mean(abs(envelope[inner] - 1.0)) < 0.02)
})

test_that(".hilbert_analytic handles edge cases", {
  # Empty input
  expect_length(PhysioCrossModal:::.hilbert_analytic(numeric(0)), 0)
  # Single element
  result <- PhysioCrossModal:::.hilbert_analytic(1.0)
  expect_length(result, 1)
  # Even length
  result_even <- PhysioCrossModal:::.hilbert_analytic(rnorm(100))
  expect_length(result_even, 100)
  # Odd length
  result_odd <- PhysioCrossModal:::.hilbert_analytic(rnorm(101))
  expect_length(result_odd, 101)
})

# ---- Bandpass Filter --------------------------------------------------------

test_that(".bandpass_filter preserves in-band signal", {
  sr <- 500
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)

  # 10 Hz signal through 5-15 Hz passband should be preserved
  x <- sin(2 * pi * 10 * t)
  filtered <- PhysioCrossModal:::.bandpass_filter(x, low = 5, high = 15, sr = sr)

  expect_length(filtered, n)

  # After transient settles (skip first/last 10%), power should be similar
  inner <- seq(floor(n * 0.1), floor(n * 0.9))
  power_original <- mean(x[inner]^2)
  power_filtered <- mean(filtered[inner]^2)
  ratio <- power_filtered / power_original

  expect_true(ratio > 0.7,
    info = paste("Power ratio:", round(ratio, 3)))
})

test_that(".bandpass_filter attenuates out-of-band signal", {
  sr <- 500
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)

  # 10 Hz signal through 20-30 Hz passband should be attenuated
  x <- sin(2 * pi * 10 * t)
  filtered <- PhysioCrossModal:::.bandpass_filter(x, low = 20, high = 30, sr = sr)

  inner <- seq(floor(n * 0.1), floor(n * 0.9))
  power_original <- mean(x[inner]^2)
  power_filtered <- mean(filtered[inner]^2)
  attenuation <- power_filtered / power_original

  expect_true(attenuation < 0.1,
    info = paste("Attenuation ratio:", round(attenuation, 5),
                 "(should be < 0.1)"))
})

test_that(".bandpass_filter validates inputs", {
  x <- rnorm(1000)
  expect_error(PhysioCrossModal:::.bandpass_filter(x, -1, 10, sr = 500),
               "positive")
  expect_error(PhysioCrossModal:::.bandpass_filter(x, 20, 10, sr = 500),
               "less than")
  expect_error(PhysioCrossModal:::.bandpass_filter(x, 10, 260, sr = 500),
               "Nyquist")
})

test_that(".bandpass_filter FIR type works", {
  sr <- 500
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.5 * sin(2 * pi * 50 * t)

  # Skip if signal package not available
  skip_if_not_installed("signal")

  filtered <- PhysioCrossModal:::.bandpass_filter(x, low = 5, high = 15,
                                                   sr = sr, type = "fir")
  expect_length(filtered, n)

  # The 50 Hz component should be attenuated
  inner <- seq(floor(n * 0.2), floor(n * 0.8))
  # Check via PSD that power at 50Hz is reduced
  psd_orig <- PhysioCrossModal:::.welch_psd(x[inner], nperseg = 256L, sr = sr)
  psd_filt <- PhysioCrossModal:::.welch_psd(filtered[inner], nperseg = 256L,
                                             sr = sr)
  # Find the 50 Hz bin
  idx_50 <- which.min(abs(psd_orig$frequencies - 50))
  ratio_50 <- psd_filt$psd[idx_50] / psd_orig$psd[idx_50]
  expect_true(ratio_50 < 0.1)
})

# ---- Signal Pair Extraction -------------------------------------------------

test_that(".extract_signal_pair with two numeric vectors returns them unchanged", {
  x <- rnorm(100)
  y <- rnorm(100)
  result <- PhysioCrossModal:::.extract_signal_pair(x, y, sr = 250)

  expect_type(result, "list")
  expect_equal(result$x, x)
  expect_equal(result$y, y)
  expect_equal(result$sr, 250)
})

test_that(".extract_signal_pair requires sr for numeric vectors", {
  expect_error(
    PhysioCrossModal:::.extract_signal_pair(rnorm(10), rnorm(10)),
    "sr.*required"
  )
})

test_that(".extract_signal_pair with two PhysioExperiment objects", {
  skip_if_not_installed("PhysioCore")

  sr_eeg <- 500
  sr_emg <- 500
  n <- 1000

  pe1 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n, ncol = 2)),
    samplingRate = sr_eeg
  )
  pe2 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 3), nrow = n, ncol = 3)),
    samplingRate = sr_emg
  )

  result <- PhysioCrossModal:::.extract_signal_pair(pe1, pe2)

  expect_type(result, "list")
  expect_equal(result$sr, sr_eeg)
  expect_length(result$x, n)
  expect_length(result$y, n)
})

test_that(".extract_signal_pair with PE objects and different sr resamples y", {
  skip_if_not_installed("PhysioCore")

  sr1 <- 500
  sr2 <- 1000
  n_sec <- 2
  n1 <- sr1 * n_sec
  n2 <- sr2 * n_sec

  pe1 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n1 * 2), nrow = n1, ncol = 2)),
    samplingRate = sr1
  )
  pe2 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n2 * 2), nrow = n2, ncol = 2)),
    samplingRate = sr2
  )

  result <- PhysioCrossModal:::.extract_signal_pair(pe1, pe2)

  expect_equal(result$sr, sr1)
  # After resampling, both should be the same length
  expect_equal(length(result$x), length(result$y))
})

test_that(".extract_signal_pair with channels_x/channels_y", {
  skip_if_not_installed("PhysioCore")

  n <- 500
  data1 <- matrix(1:4, nrow = n, ncol = 2, byrow = TRUE)
  data1[, 1] <- seq_len(n)
  data1[, 2] <- seq_len(n) * 10

  pe1 <- PhysioExperiment(
    assays = list(raw = data1),
    samplingRate = 250
  )

  # Extract channel 2 from pe1 and channel 1 from pe1
  result <- PhysioCrossModal:::.extract_signal_pair(
    pe1, pe1, channels_x = 2L, channels_y = 1L
  )

  expect_equal(result$x, as.numeric(data1[, 2]))
  expect_equal(result$y, as.numeric(data1[, 1]))
})

test_that(".extract_signal_pair rejects unsupported inputs", {
  expect_error(
    PhysioCrossModal:::.extract_signal_pair("not_a_signal"),
    "Unsupported"
  )
})

# ---- Resample Linear -------------------------------------------------------

test_that(".resample_linear resamples correctly", {
  # Create a known signal and resample
  sr_from <- 100
  sr_to <- 200
  n <- 100
  t <- seq(0, (n - 1) / sr_from, length.out = n)
  x <- sin(2 * pi * 5 * t)  # 5 Hz sine

  resampled <- PhysioCrossModal:::.resample_linear(x, sr_from, sr_to)

  # Should have approximately 2x the number of samples
  expected_n <- as.integer(round((n - 1) / sr_from * sr_to)) + 1L
  expect_equal(length(resampled), expected_n)

  # The resampled signal should still look like a 5Hz sine.
  # Linear interpolation of a sinusoid introduces some error at steep parts,
  # so we use a relaxed tolerance.
  t2 <- seq(0, (n - 1) / sr_from, length.out = expected_n)
  expected <- sin(2 * pi * 5 * t2)
  expect_equal(resampled, expected, tolerance = 0.02)
})

# ---- Integration tests (combining utilities) --------------------------------

test_that("Welch PSD and CSD are consistent for pure sinusoid", {
  sr <- 500
  n <- 4000
  t <- seq(0, (n - 1) / sr, length.out = n)
  freq <- 25
  x <- sin(2 * pi * freq * t)

  psd <- PhysioCrossModal:::.welch_psd(x, nperseg = 512L, sr = sr)
  csd <- PhysioCrossModal:::.welch_csd(x, x, nperseg = 512L, sr = sr)

  # Both should show peak at 25 Hz
  psd_peak <- psd$frequencies[which.max(psd$psd)]
  csd_peak <- csd$frequencies[which.max(Mod(csd$csd))]
  expect_equal(psd_peak, csd_peak, tolerance = 2)
  expect_true(abs(psd_peak - freq) < 2)
})

test_that("Filter + Hilbert pipeline works for amplitude modulated signal", {
  sr <- 500
  n <- 5000
  t <- seq(0, (n - 1) / sr, length.out = n)

  # AM signal: carrier at 40 Hz, modulation at 5 Hz
  carrier <- sin(2 * pi * 40 * t)
  modulation <- 0.5 + 0.5 * cos(2 * pi * 5 * t)
  x <- modulation * carrier

  # Bandpass around carrier
  filtered <- PhysioCrossModal:::.bandpass_filter(x, low = 30, high = 50, sr = sr)

  # Extract envelope
  envelope <- PhysioCrossModal:::.hilbert_envelope(filtered)

  # Demean the envelope to remove the DC component before spectral analysis
  envelope_dm <- envelope - mean(envelope)

  # Envelope should modulate at ~5 Hz. Use larger nperseg for better resolution.
  env_psd <- PhysioCrossModal:::.welch_psd(envelope_dm, nperseg = 512L, sr = sr)

  # Find peak frequency (DC component is already removed by demeaning,
  # but skip the DC bin just in case)
  search_idx <- which(env_psd$frequencies >= 2 & env_psd$frequencies <= 20)
  peak_idx <- search_idx[which.max(env_psd$psd[search_idx])]
  peak_freq <- env_psd$frequencies[peak_idx]
  expect_true(abs(peak_freq - 5) < 3,
    info = paste("Envelope modulation peak at", peak_freq, "Hz, expected ~5 Hz"))
})

# ---- .lowpass_filter tests ---------------------------------------------------

test_that(".lowpass_filter preserves low-frequency content", {
  sr <- 500
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- sin(2 * pi * 5 * t)  # 5 Hz signal

  filtered <- PhysioCrossModal:::.lowpass_filter(x, cutoff = 50, sr = sr)

  inner <- seq(floor(n * 0.1), floor(n * 0.9))
  power_ratio <- mean(filtered[inner]^2) / mean(x[inner]^2)
  expect_true(power_ratio > 0.8,
    info = paste("Power ratio:", round(power_ratio, 3)))
})

test_that(".lowpass_filter attenuates high-frequency content", {
  sr <- 500
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- sin(2 * pi * 100 * t)  # 100 Hz signal

  filtered <- PhysioCrossModal:::.lowpass_filter(x, cutoff = 30, sr = sr)

  inner <- seq(floor(n * 0.1), floor(n * 0.9))
  power_ratio <- mean(filtered[inner]^2) / mean(x[inner]^2)
  expect_true(power_ratio < 0.1,
    info = paste("Power ratio:", round(power_ratio, 5)))
})

test_that(".lowpass_filter returns input unchanged when cutoff >= Nyquist", {
  sr <- 500
  x <- rnorm(1000)
  result <- PhysioCrossModal:::.lowpass_filter(x, cutoff = 250, sr = sr)
  expect_equal(result, x)
})

# ---- .resolve_window tests ---------------------------------------------------

test_that(".resolve_window passes numeric vector through", {
  win <- c(0.5, 1.0, 0.5)
  result <- PhysioCrossModal:::.resolve_window(win, 3)
  expect_equal(result, win)
})

test_that(".resolve_window resolves 'hann' alias", {
  result <- PhysioCrossModal:::.resolve_window("hann", 4)
  expected <- PhysioCrossModal:::.window_hanning(4)
  expect_equal(result, expected)
})

test_that(".resolve_window errors on unknown window name", {
  expect_error(
    PhysioCrossModal:::.resolve_window("kaiser", 10),
    "Unknown window type"
  )
})

test_that(".resolve_window errors on numeric window with wrong length", {
  expect_error(
    PhysioCrossModal:::.resolve_window(c(1, 2, 3), 5),
    "length equal to nperseg"
  )
})

# ---- .bandpass_filter near-boundary clamping warning -------------------------

test_that(".bandpass_filter warns when cutoff near boundary requires clamping", {
  sr <- 500
  x <- rnorm(2000)
  # low = 0.5 Hz with Nyquist = 250, normalized = 0.002 which clamps to 0.01
  expect_warning(
    PhysioCrossModal:::.bandpass_filter(x, low = 0.5, high = 100, sr = sr),
    "clamped"
  )
})

# ---- Numerical validation: Parseval's theorem --------------------------------

test_that("Welch PSD satisfies Parseval's theorem (sum(psd)*df ≈ var(x))", {
  set.seed(42)
  sr <- 1000
  n <- 4000
  x <- rnorm(n)

  result <- PhysioCrossModal:::.welch_psd(x, nperseg = 512L, sr = sr)
  df <- result$frequencies[2] - result$frequencies[1]
  psd_power <- sum(result$psd) * df
  signal_var <- stats::var(x)

  # Parseval: integral of PSD ≈ variance
  ratio <- psd_power / signal_var
  expect_true(abs(ratio - 1) < 0.3,
    info = sprintf("Parseval ratio = %.3f (expect ~1)", ratio))
})

# ---- Self-coherence ----------------------------------------------------------

test_that("coherence(x, x) = 1.0 at all frequencies", {
  sr <- 500
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.5 * sin(2 * pi * 30 * t) + 0.2 * rnorm(n)

  result <- coherence(x, x, sr = sr, nperseg = 256L)

  # Self-coherence should be 1.0 everywhere
  expect_true(all(result$coherence > 0.99),
    info = paste("Min self-coherence:", min(result$coherence)))
})

# ---- CSD symmetry ------------------------------------------------------------

test_that("CSD(x,y) = conj(CSD(y,x)) (Hermitian symmetry)", {
  set.seed(42)
  sr <- 500
  n <- 2000
  x <- rnorm(n)
  y <- rnorm(n)

  csd_xy <- PhysioCrossModal:::.welch_csd(x, y, nperseg = 256L, sr = sr)
  csd_yx <- PhysioCrossModal:::.welch_csd(y, x, nperseg = 256L, sr = sr)

  expect_equal(csd_xy$csd, Conj(csd_yx$csd), tolerance = 1e-10)
})

# ---- White noise PSD approximately flat --------------------------------------

test_that("PSD of white noise is approximately flat", {
  set.seed(42)
  sr <- 1000
  n <- 10000
  x <- rnorm(n)

  result <- PhysioCrossModal:::.welch_psd(x, nperseg = 512L, sr = sr)

  # Skip DC and Nyquist bins
  inner <- result$psd[2:(length(result$psd) - 1)]
  cv <- stats::sd(inner) / mean(inner)  # coefficient of variation
  expect_true(cv < 0.5,
    info = paste("CV of white noise PSD:", round(cv, 3)))
})

# ---- Edge cases: .welch_psd -------------------------------------------------

test_that(".welch_psd with n=2 (minimum) works", {
  x <- c(1, -1)
  result <- PhysioCrossModal:::.welch_psd(x, nperseg = 256L, sr = 100)
  expect_type(result, "list")
  expect_true(length(result$psd) > 0)
})

test_that(".welch_psd with nperseg >= signal length adapts", {
  x <- rnorm(100)
  result <- PhysioCrossModal:::.welch_psd(x, nperseg = 200L, sr = 100)
  expect_type(result, "list")
  expect_true(length(result$psd) > 0)
})

# ---- Edge case: .hilbert_analytic with all-zeros -----------------------------

test_that(".hilbert_analytic with all-zeros returns all-zeros", {
  x <- rep(0, 100)
  result <- PhysioCrossModal:::.hilbert_analytic(x)
  expect_length(result, 100)
  expect_true(all(Mod(result) < .Machine$double.eps))
})

# ---- Edge case: .bandpass_filter with short signal ---------------------------

test_that(".bandpass_filter handles signal shorter than filter order", {
  sr <- 500
  x <- sin(2 * pi * 10 * seq(0, 0.02, length.out = 10))

  # Should fall back to FFT method without error
  result <- PhysioCrossModal:::.bandpass_filter(x, low = 5, high = 100, sr = sr)
  expect_length(result, length(x))
  expect_true(is.numeric(result))
})
