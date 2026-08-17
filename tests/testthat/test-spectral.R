# Tests for spectral coupling functions (coupling-spectral.R)

# ---- coherence() -------------------------------------------------------------

test_that("coherence detects peak at 20Hz for coupled sinusoids", {
  sigs <- make_coupled_signals(
    sr1 = 500, sr2 = 500, n_sec = 10,
    coupling_freq = 20, coupling_strength = 0.8, noise_level = 0.2
  )

  result <- coherence(sigs$x, sigs$y, sr = sigs$sr_x, nperseg = 512L)

  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true("confidence_limit" %in% names(result))
  expect_true("n_segments" %in% names(result))

  # Find the coherence at the coupling frequency (20 Hz)
  idx_20 <- which.min(abs(result$frequencies - 20))
  expect_true(result$coherence[idx_20] > 0.5,
    info = paste("Coherence at 20 Hz:", result$coherence[idx_20],
                 "-- expected > 0.5"))
})

test_that("coherence is low for independent noise signals", {
  set.seed(123)
  n <- 5000
  x <- rnorm(n)
  y <- rnorm(n)

  result <- coherence(x, y, sr = 500, nperseg = 256L)

  # For independent signals, coherence should be very low everywhere
  expect_true(max(result$coherence) < 0.4,
    info = paste("Max coherence for independent noise:",
                 max(result$coherence), "-- expected < 0.4"))
  # Mean coherence should be even lower
  expect_true(mean(result$coherence) < 0.2,
    info = paste("Mean coherence for independent noise:",
                 mean(result$coherence)))
})

test_that("coherence freq_range filtering works", {
  sigs <- make_coupled_signals(
    sr1 = 500, sr2 = 500, n_sec = 10,
    coupling_freq = 20, coupling_strength = 0.8, noise_level = 0.2
  )

  result_full <- coherence(sigs$x, sigs$y, sr = sigs$sr_x, nperseg = 256L)
  result_filt <- coherence(sigs$x, sigs$y, sr = sigs$sr_x, nperseg = 256L,
                           freq_range = c(10, 30))

  # Filtered result should only have frequencies in [10, 30]
  expect_true(all(result_filt$frequencies >= 10))
  expect_true(all(result_filt$frequencies <= 30))

  # Filtered result should have fewer frequency bins than full result
  expect_true(length(result_filt$frequencies) < length(result_full$frequencies))

  # Coherence values in the filtered range should match
  idx_full <- which(result_full$frequencies >= 10 & result_full$frequencies <= 30)
  expect_equal(result_filt$coherence, result_full$coherence[idx_full])
  expect_equal(result_filt$frequencies, result_full$frequencies[idx_full])
})

test_that("coherence returns positive confidence_limit", {
  set.seed(42)
  x <- rnorm(2000)
  y <- rnorm(2000)

  result <- coherence(x, y, sr = 500, nperseg = 256L)

  expect_true(is.numeric(result$confidence_limit))
  expect_length(result$confidence_limit, 1)
  expect_true(result$confidence_limit > 0,
    info = paste("Confidence limit:", result$confidence_limit))
  expect_true(result$confidence_limit < 1)
})

test_that("coherence n_segments is returned correctly", {
  set.seed(42)
  n <- 2000
  x <- rnorm(n)
  y <- rnorm(n)

  result <- coherence(x, y, sr = 500, nperseg = 256L)

  expect_true(is.integer(result$n_segments) || is.numeric(result$n_segments))
  expect_true(result$n_segments >= 1)

  # With nperseg=256 and noverlap=128, step=128, segments ~ (2000-128)/128 ~ 14
  expect_true(result$n_segments > 10)
})

test_that("coherence works with two numeric vectors", {
  set.seed(42)
  sr <- 250
  t <- seq(0, 4, length.out = sr * 4)
  x <- sin(2 * pi * 15 * t) + 0.3 * rnorm(length(t))
  y <- 0.7 * sin(2 * pi * 15 * t) + 0.3 * rnorm(length(t))

  result <- coherence(x, y, sr = sr, nperseg = 128L)

  expect_type(result, "list")
  expect_length(result$coherence, length(result$frequencies))

  # Should detect coupling at 15 Hz

  idx_15 <- which.min(abs(result$frequencies - 15))
  expect_true(result$coherence[idx_15] > 0.3)
})

test_that("coherence works with two PhysioExperiment objects", {
  skip_if_not_installed("PhysioCore")

  pair <- make_eeg_emg(n_sec = 10, sr_eeg = 500, sr_emg = 500,
                       cmc_freq = 20, cmc_strength = 0.8)

  result <- coherence(pair$eeg, pair$emg, nperseg = 512L,
                      channels_x = 1L, channels_y = 1L)

  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true("confidence_limit" %in% names(result))

  # Should detect coupling at 20 Hz in the first channels
  idx_20 <- which.min(abs(result$frequencies - 20))
  expect_true(result$coherence[idx_20] > 0.3,
    info = paste("Coherence at 20 Hz (PE input):", result$coherence[idx_20]))
})

test_that("coherence values are bounded in [0, 1]", {
  sigs <- make_coupled_signals(sr1 = 500, sr2 = 500, n_sec = 5)

  result <- coherence(sigs$x, sigs$y, sr = sigs$sr_x, nperseg = 256L)

  expect_true(all(result$coherence >= 0))
  expect_true(all(result$coherence <= 1))
})

test_that("coherence of identical signals is near 1.0", {
  sr <- 500
  n <- 5000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.1 * rnorm(n)

  result <- coherence(x, x, sr = sr, nperseg = 256L)

  # Coherence of a signal with itself should be very high everywhere
  # (not exactly 1.0 because Welch averaging introduces some variability)
  expect_true(mean(result$coherence) > 0.95,
    info = paste("Mean self-coherence:", mean(result$coherence)))
})

test_that("coherence rejects multitaper method (not yet implemented)", {
  expect_error(
    coherence(rnorm(1000), rnorm(1000), sr = 500, method = "multitaper"),
    "not yet implemented"
  )
})

# ---- crossSpectrum() ---------------------------------------------------------

test_that("crossSpectrum returns complex CSD with correct frequency axis", {
  set.seed(42)
  sr <- 500
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.1 * rnorm(n)
  y <- sin(2 * pi * 10 * t + pi / 4) + 0.1 * rnorm(n)

  result <- crossSpectrum(x, y, sr = sr, nperseg = 256L)

  expect_type(result, "list")
  expect_true("csd" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true("magnitude" %in% names(result))
  expect_true("phase" %in% names(result))

  # CSD should be complex
  expect_true(is.complex(result$csd))

  # Frequency axis should span [0, Nyquist]
  expect_equal(result$frequencies[1], 0)
  expect_equal(result$frequencies[length(result$frequencies)], sr / 2)

  # All vectors should have the same length
  n_freqs <- length(result$frequencies)
  expect_length(result$csd, n_freqs)
  expect_length(result$magnitude, n_freqs)
  expect_length(result$phase, n_freqs)
})

test_that("crossSpectrum magnitude and phase are correct", {
  set.seed(42)
  sr <- 500
  n <- 4000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- sin(2 * pi * 25 * t) + 0.05 * rnorm(n)
  y <- sin(2 * pi * 25 * t) + 0.05 * rnorm(n)

  result <- crossSpectrum(x, y, sr = sr, nperseg = 512L)

  # Magnitude should equal Mod(csd)
  expect_equal(result$magnitude, Mod(result$csd))

  # Phase should equal Arg(csd)
  expect_equal(result$phase, Arg(result$csd))

  # Magnitude should peak near 25 Hz
  peak_idx <- which.max(result$magnitude)
  peak_freq <- result$frequencies[peak_idx]
  expect_true(abs(peak_freq - 25) < 2,
    info = paste("CSD magnitude peak at", peak_freq, "Hz, expected ~25 Hz"))

  # For two nearly identical signals, phase at peak should be near 0
  expect_true(abs(result$phase[peak_idx]) < 0.3,
    info = paste("Phase at peak:", result$phase[peak_idx], "expected ~0"))
})

test_that("crossSpectrum detects phase difference", {
  sr <- 500
  n <- 5000
  t <- seq(0, (n - 1) / sr, length.out = n)
  phase_shift <- pi / 2  # 90 degrees

  x <- sin(2 * pi * 10 * t)
  y <- sin(2 * pi * 10 * t + phase_shift)

  result <- crossSpectrum(x, y, sr = sr, nperseg = 512L)

  # Find the 10 Hz bin
  idx_10 <- which.min(abs(result$frequencies - 10))

  # Phase at the 10 Hz bin should be approximately pi/2
  # Note: CSD phase convention is Arg(X * conj(Y)) = phase_x - phase_y
  # So for y shifted by +pi/2, phase should be approximately -pi/2
  expect_true(abs(abs(result$phase[idx_10]) - pi / 2) < 0.3,
    info = paste("Phase at 10 Hz:", result$phase[idx_10],
                 "expected near +/-pi/2"))
})

test_that("crossSpectrum works with numeric vectors", {
  set.seed(42)
  x <- rnorm(1000)
  y <- rnorm(1000)

  result <- crossSpectrum(x, y, sr = 250, nperseg = 128L)

  expect_type(result, "list")
  expect_true(is.complex(result$csd))
  expect_true(is.numeric(result$magnitude))
  expect_true(is.numeric(result$phase))
})

test_that("crossSpectrum works with PhysioExperiment objects", {
  skip_if_not_installed("PhysioCore")

  pair <- make_eeg_emg(n_sec = 5, sr_eeg = 500, sr_emg = 500,
                       cmc_freq = 20, cmc_strength = 0.7)

  result <- crossSpectrum(pair$eeg, pair$emg,
                          channels_x = 1L, channels_y = 1L,
                          nperseg = 256L)

  expect_type(result, "list")
  expect_true(is.complex(result$csd))

  # Should show magnitude peak near 20 Hz
  peak_idx <- which.max(result$magnitude)
  peak_freq <- result$frequencies[peak_idx]
  expect_true(abs(peak_freq - 20) < 5,
    info = paste("CSD peak at", peak_freq, "Hz for PE input, expected ~20 Hz"))
})

test_that("crossSpectrum of signal with itself matches PSD shape", {
  sr <- 500
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)
  x <- sin(2 * pi * 15 * t) + 0.1 * rnorm(n)

  cs_result <- crossSpectrum(x, x, sr = sr, nperseg = 256L)

  # CSD of x with itself should be purely real (imaginary part ~ 0)
  imag_ratio <- max(abs(Im(cs_result$csd))) / max(abs(Re(cs_result$csd)))
  expect_true(imag_ratio < 1e-10,
    info = paste("Imaginary/Real ratio for self-CSD:", imag_ratio))

  # Magnitude should peak at 15 Hz
  peak_idx <- which.max(cs_result$magnitude)
  peak_freq <- cs_result$frequencies[peak_idx]
  expect_true(abs(peak_freq - 15) < 2)
})
