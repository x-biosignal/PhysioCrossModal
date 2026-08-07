# Tests for multitaper spectral estimation

# ---- .dpss_tapers ------------------------------------------------------------

test_that(".dpss_tapers returns correct dimensions", {
  tapers <- PhysioCrossModal:::.dpss_tapers(100, nw = 4, k = 7)
  expect_true(is.matrix(tapers))
  expect_equal(dim(tapers), c(100, 7))
})

test_that(".dpss_tapers default k = 2*nw - 1", {
  tapers <- PhysioCrossModal:::.dpss_tapers(64, nw = 3)
  expect_equal(ncol(tapers), 5)
})

test_that(".dpss_tapers are approximately orthonormal", {
  tapers <- PhysioCrossModal:::.dpss_tapers(100, nw = 4, k = 7)
  # Check orthogonality: t(tapers) %*% tapers should be close to identity
  gram <- crossprod(tapers)
  expect_equal(gram, diag(7), tolerance = 0.1)
})

# ---- .multitaper_psd ---------------------------------------------------------

test_that(".multitaper_psd peak at known frequency", {
  sr <- 256
  n <- sr * 4
  t_vec <- seq(0, 4 - 1/sr, length.out = n)
  x <- sin(2 * pi * 20 * t_vec)

  result <- PhysioCrossModal:::.multitaper_psd(x, sr, nw = 4)

  peak_idx <- which.max(result$psd)
  peak_freq <- result$frequencies[peak_idx]
  expect_true(abs(peak_freq - 20) < 2,
    info = paste("Peak at", peak_freq, "Hz, expected ~20 Hz"))
})

test_that(".multitaper_psd returns correct structure", {
  x <- rnorm(128)
  result <- PhysioCrossModal:::.multitaper_psd(x, sr = 100, nw = 4)

  expect_type(result, "list")
  expect_true("psd" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_length(result$psd, length(result$frequencies))
  expect_true(all(result$psd >= 0))
})

# ---- .multitaper_csd ---------------------------------------------------------

test_that(".multitaper_csd returns correct structure", {
  x <- rnorm(200)
  y <- rnorm(200)
  result <- PhysioCrossModal:::.multitaper_csd(x, y, sr = 100, nw = 4)

  expect_type(result, "list")
  expect_true("csd" %in% names(result))
  expect_true("psd_x" %in% names(result))
  expect_true("psd_y" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true(is.complex(result$csd))
})

test_that(".multitaper_csd self CSD equals PSD", {
  set.seed(42)
  x <- rnorm(200)
  result <- PhysioCrossModal:::.multitaper_csd(x, x, sr = 100, nw = 4)

  # CSD(x,x) should be real and equal to PSD(x)
  expect_true(all(abs(Im(result$csd)) < max(Mod(result$csd)) * 0.01))
  expect_equal(Re(result$csd), result$psd_x, tolerance = 1e-10)
})

# ---- multitaperCoherence -----------------------------------------------------

test_that("multitaperCoherence self-coherence is approximately 1", {
  set.seed(42)
  sr <- 256
  n <- sr * 4
  t_vec <- seq(0, 4 - 1/sr, length.out = n)
  x <- sin(2 * pi * 20 * t_vec) + 0.1 * rnorm(n)

  result <- multitaperCoherence(x, x, sr = sr, nw = 4)

  # Self-coherence should be very close to 1
  expect_true(mean(result$coherence) > 0.95)
})

test_that("multitaperCoherence peak at coupling frequency", {
  set.seed(42)
  sr <- 256
  n <- sr * 4
  t_vec <- seq(0, 4 - 1/sr, length.out = n)
  x <- sin(2 * pi * 20 * t_vec) + 0.2 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 20 * t_vec) + 0.2 * rnorm(n)

  result <- multitaperCoherence(x, y, sr = sr, nw = 4)

  peak_idx <- which.max(result$coherence)
  peak_freq <- result$frequencies[peak_idx]
  expect_true(abs(peak_freq - 20) < 3,
    info = paste("Peak at", peak_freq, "Hz, expected ~20 Hz"))
})

test_that("multitaperCoherence independent signals below threshold", {
  set.seed(42)
  x <- rnorm(1024)
  y <- rnorm(1024)

  result <- multitaperCoherence(x, y, sr = 256, nw = 4)

  mean_coh <- mean(result$coherence)
  expect_true(mean_coh < 0.5,
    info = paste("Mean coherence =", round(mean_coh, 3)))
})

test_that("multitaperCoherence values in [0, 1]", {
  set.seed(42)
  x <- rnorm(512)
  y <- rnorm(512)

  result <- multitaperCoherence(x, y, sr = 256, nw = 4)

  expect_true(all(result$coherence >= 0 & result$coherence <= 1))
})

test_that("multitaperCoherence with PE input", {
  skip_if_not_installed("PhysioCore")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  t_vec <- seq(0, 2, length.out = n)
  data1 <- matrix(sin(2 * pi * 10 * t_vec) + 0.2 * rnorm(n), ncol = 1)
  data2 <- matrix(0.8 * sin(2 * pi * 10 * t_vec) + 0.2 * rnorm(n), ncol = 1)

  pe1 <- PhysioExperiment(assays = list(raw = data1), samplingRate = sr)
  pe2 <- PhysioExperiment(assays = list(raw = data2), samplingRate = sr)

  result <- multitaperCoherence(pe1, pe2)
  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
})

test_that("multitaperCoherence with MPE input", {
  skip_if_not_installed("PhysioCore")

  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- multitaperCoherence(mpe, modality_x = "EEG", modality_y = "EMG")
  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
})

test_that("multitaperCoherence freq_range filtering", {
  set.seed(42)
  x <- rnorm(512)
  y <- rnorm(512)

  result <- multitaperCoherence(x, y, sr = 256, freq_range = c(10, 50))

  expect_true(all(result$frequencies >= 10 & result$frequencies <= 50))
})

test_that("multitaperCoherence confidence_limit is reasonable", {
  set.seed(42)
  x <- rnorm(512)
  y <- rnorm(512)

  result <- multitaperCoherence(x, y, sr = 256, nw = 4)

  expect_true(result$confidence_limit > 0 && result$confidence_limit < 1)
})

test_that("couplingAnalysis dispatches multitaper_coherence", {
  set.seed(42)
  sr <- 256
  n <- sr * 2
  t_vec <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 20 * t_vec) + 0.2 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 20 * t_vec) + 0.2 * rnorm(n)

  result <- couplingAnalysis(x, y, method = "multitaper_coherence", sr = sr)

  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
  expect_true("frequencies" %in% names(result))
})
