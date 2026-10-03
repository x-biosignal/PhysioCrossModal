# Tests for wavelet-based coupling (coupling-wavelet.R)

# ---- .morlet_wavelet_transform -----------------------------------------------

test_that(".morlet_wavelet_transform returns correct dimensions", {
  x <- rnorm(500)
  freqs <- c(5, 10, 20)
  result <- PhysioCrossModal:::.morlet_wavelet_transform(x, freqs, n_cycles = 7,
                                                          sr = 250)
  expect_true(is.matrix(result))
  expect_true(is.complex(result))
  expect_equal(dim(result), c(500, 3))
})

test_that(".morlet_wavelet_transform shows peak at signal frequency", {
  sr <- 200
  n <- sr * 2
  t <- seq(0, 2 - 1/sr, length.out = n)
  x <- sin(2 * pi * 10 * t)

  freqs <- seq(5, 20, by = 1)
  W <- PhysioCrossModal:::.morlet_wavelet_transform(x, freqs, n_cycles = 7,
                                                     sr = sr)
  # Average power across time (skip edges)
  inner <- seq(floor(n * 0.1), floor(n * 0.9))
  avg_power <- colMeans(Mod(W[inner, ])^2)

  # Peak should be at 10 Hz
  peak_freq <- freqs[which.max(avg_power)]
  expect_true(abs(peak_freq - 10) <= 2,
    info = paste("Peak at", peak_freq, "Hz, expected ~10 Hz"))
})

# ---- waveletCoherence --------------------------------------------------------

test_that("waveletCoherence returns correct structure", {
  set.seed(42)
  sr <- 200
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t) + 0.2 * rnorm(n)

  result <- waveletCoherence(x, y, sr = sr, frequencies = seq(5, 20))

  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
  expect_true("phase" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true("times" %in% names(result))

  expect_equal(ncol(result$coherence), length(result$frequencies))
  expect_equal(nrow(result$coherence), length(result$times))
  expect_equal(dim(result$phase), dim(result$coherence))
})

test_that("waveletCoherence values are in [0, 1]", {
  set.seed(42)
  sr <- 200
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t) + 0.2 * rnorm(n)

  result <- waveletCoherence(x, y, sr = sr, frequencies = seq(5, 20))

  expect_true(all(result$coherence >= 0 & result$coherence <= 1))
})

test_that("waveletCoherence peaks at coupling frequency", {
  set.seed(42)
  sr <- 200
  n <- sr * 4
  t <- seq(0, 4, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t) + 0.2 * rnorm(n)

  freqs <- seq(5, 20, by = 1)
  result <- waveletCoherence(x, y, sr = sr, frequencies = freqs)

  # Average coherence over time (skip edges)
  inner <- seq(floor(n * 0.2), floor(n * 0.8))
  avg_coh <- colMeans(result$coherence[inner, ])

  peak_freq <- freqs[which.max(avg_coh)]
  expect_true(abs(peak_freq - 10) <= 3,
    info = paste("Peak coherence at", peak_freq, "Hz, expected ~10 Hz"))
})

test_that("waveletCoherence low for independent signals", {
  set.seed(42)
  sr <- 200
  n <- sr * 2
  x <- rnorm(n)
  y <- rnorm(n)

  result <- waveletCoherence(x, y, sr = sr, frequencies = seq(5, 20))

  # Mean coherence should be relatively low
  mean_coh <- mean(result$coherence)
  expect_true(mean_coh < 0.6,
    info = paste("Mean coherence =", round(mean_coh, 3),
                 "for independent signals"))
})

test_that("waveletCoherence works with PhysioExperiment input", {
  skip_if_not_installed("PhysioCore")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  data1 <- matrix(sin(2 * pi * 10 * t) + 0.2 * rnorm(n), ncol = 1)
  data2 <- matrix(0.8 * sin(2 * pi * 10 * t) + 0.2 * rnorm(n), ncol = 1)

  pe1 <- PhysioExperiment(assays = list(raw = data1), samplingRate = sr)
  pe2 <- PhysioExperiment(assays = list(raw = data2), samplingRate = sr)

  result <- waveletCoherence(pe1, pe2, frequencies = seq(5, 15))

  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
})

test_that("waveletCoherence works with MultiPhysioExperiment input", {
  skip_if_not_installed("PhysioCore")

  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- waveletCoherence(mpe, modality_x = "EEG", modality_y = "EMG",
                             frequencies = seq(10, 30))

  expect_type(result, "list")
  expect_true(ncol(result$coherence) == length(seq(10, 30)))
})

test_that("waveletCoherence detects time-varying coupling", {
  set.seed(42)
  sr <- 200
  n <- sr * 4
  t <- seq(0, 4, length.out = n)

  # Coupling only in the first half
  coupling <- c(rep(1, n %/% 2), rep(0, n - n %/% 2))
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- coupling * sin(2 * pi * 10 * t) + 0.3 * rnorm(n)

  result <- waveletCoherence(x, y, sr = sr, frequencies = c(10))

  # Coherence should be higher in the first half than the second
  mid <- n %/% 2
  coh_first <- mean(result$coherence[floor(mid * 0.3):floor(mid * 0.7), 1])
  coh_second <- mean(result$coherence[floor(mid * 1.3):floor(mid * 1.7), 1])

  expect_true(coh_first > coh_second,
    info = sprintf("First half coh = %.3f, second half coh = %.3f",
                   coh_first, coh_second))
})

# ---- waveletPLV --------------------------------------------------------------

test_that("waveletPLV returns correct structure", {
  set.seed(42)
  sr <- 200
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- sin(2 * pi * 10 * t + pi / 4) + 0.2 * rnorm(n)

  result <- waveletPLV(x, y, sr = sr, frequencies = seq(5, 20))

  expect_type(result, "list")
  expect_true("plv" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true("times" %in% names(result))

  expect_equal(ncol(result$plv), length(result$frequencies))
  expect_equal(nrow(result$plv), length(result$times))
})

test_that("waveletPLV values are in [0, 1]", {
  set.seed(42)
  sr <- 200
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- sin(2 * pi * 10 * t + pi / 4) + 0.2 * rnorm(n)

  result <- waveletPLV(x, y, sr = sr, frequencies = seq(5, 20))

  expect_true(all(result$plv >= 0 & result$plv <= 1))
})

test_that("waveletPLV peaks at coupled frequency", {
  set.seed(42)
  sr <- 200
  n <- sr * 4
  t <- seq(0, 4, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- sin(2 * pi * 10 * t + pi / 4) + 0.2 * rnorm(n)

  freqs <- seq(5, 20, by = 1)
  result <- waveletPLV(x, y, sr = sr, frequencies = freqs)

  # Average PLV over time (skip edges)
  inner <- seq(floor(n * 0.2), floor(n * 0.8))
  avg_plv <- colMeans(result$plv[inner, ])

  peak_freq <- freqs[which.max(avg_plv)]
  expect_true(abs(peak_freq - 10) <= 3,
    info = paste("Peak PLV at", peak_freq, "Hz, expected ~10 Hz"))
})

test_that("waveletPLV low for independent signals", {
  set.seed(42)
  sr <- 200
  n <- sr * 2
  x <- rnorm(n)
  y <- rnorm(n)

  result <- waveletPLV(x, y, sr = sr, frequencies = seq(5, 20))

  mean_plv <- mean(result$plv)
  expect_true(mean_plv < 0.6,
    info = paste("Mean PLV =", round(mean_plv, 3),
                 "for independent signals"))
})

test_that("waveletPLV works with PhysioExperiment input", {
  skip_if_not_installed("PhysioCore")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  data1 <- matrix(sin(2 * pi * 10 * t), ncol = 1)
  data2 <- matrix(sin(2 * pi * 10 * t + pi / 4), ncol = 1)

  pe1 <- PhysioExperiment(assays = list(raw = data1), samplingRate = sr)
  pe2 <- PhysioExperiment(assays = list(raw = data2), samplingRate = sr)

  result <- waveletPLV(pe1, pe2, frequencies = seq(5, 15))

  expect_type(result, "list")
  expect_true("plv" %in% names(result))
})

# ---- .compute_coi -------------------------------------------------------------

test_that(".compute_coi returns numeric vector of correct length", {
  coi <- PhysioCrossModal:::.compute_coi(500, sr = 250, n_cycles = 7)
  expect_length(coi, 500)
  expect_true(is.numeric(coi))
})

test_that(".compute_coi is symmetric", {
  coi <- PhysioCrossModal:::.compute_coi(200, sr = 100, n_cycles = 7)
  # Start and end values should be equal (symmetric about center)
  expect_equal(coi[1], coi[200], tolerance = 1e-10)
  expect_equal(coi[10], coi[191], tolerance = 1e-10)
})

test_that(".compute_coi is highest at edges, lowest in center", {
  coi <- PhysioCrossModal:::.compute_coi(200, sr = 100, n_cycles = 7)
  center <- coi[100]
  edge <- coi[1]
  expect_true(edge > center)
})

test_that(".compute_coi values are positive", {
  coi <- PhysioCrossModal:::.compute_coi(100, sr = 50, n_cycles = 7)
  expect_true(all(coi > 0))
})

test_that("waveletCoherence return includes coi", {
  set.seed(42)
  sr <- 200
  n <- sr * 2
  t_vec <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t_vec) + 0.2 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t_vec) + 0.2 * rnorm(n)

  result <- waveletCoherence(x, y, sr = sr, frequencies = seq(5, 20))
  expect_true("coi" %in% names(result))
  expect_length(result$coi, n)
  expect_true(all(result$coi > 0))
})

test_that("waveletPLV return includes coi", {
  set.seed(42)
  sr <- 200
  n <- sr * 2
  t_vec <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t_vec) + 0.2 * rnorm(n)
  y <- sin(2 * pi * 10 * t_vec + pi / 4) + 0.2 * rnorm(n)

  result <- waveletPLV(x, y, sr = sr, frequencies = seq(5, 20))
  expect_true("coi" %in% names(result))
  expect_length(result$coi, n)
  expect_true(all(result$coi > 0))
})

# ---- plotWaveletCoherence ----------------------------------------------------

test_that("plotWaveletCoherence returns ggplot object", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  t_vec <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t_vec) + 0.3 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t_vec) + 0.3 * rnorm(n)

  result <- waveletCoherence(x, y, sr = sr, frequencies = seq(5, 20))
  p <- plotWaveletCoherence(result)
  expect_s3_class(p, "ggplot")
})

test_that("plotWaveletCoherence works with waveletPLV result", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  t_vec <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t_vec) + 0.3 * rnorm(n)
  y <- sin(2 * pi * 10 * t_vec + pi / 4) + 0.3 * rnorm(n)

  result <- waveletPLV(x, y, sr = sr, frequencies = seq(5, 20))
  p <- plotWaveletCoherence(result)
  expect_s3_class(p, "ggplot")
})

test_that("plotWaveletCoherence COI overlay renders without error", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  t_vec <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t_vec) + 0.3 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t_vec) + 0.3 * rnorm(n)

  result <- waveletCoherence(x, y, sr = sr, frequencies = seq(5, 20))
  expect_no_error(plotWaveletCoherence(result, show_coi = TRUE))
})

test_that("plotWaveletCoherence works without COI", {
  skip_if_not_installed("ggplot2")

  # Simulate a v0.2.0 result without coi field
  result <- list(
    coherence = matrix(runif(100), nrow = 10, ncol = 10),
    frequencies = seq(5, 14),
    times = seq(0, 0.9, by = 0.1)
  )
  p <- plotWaveletCoherence(result, show_coi = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plotWaveletCoherence custom title and fill_label", {
  skip_if_not_installed("ggplot2")

  result <- list(
    coherence = matrix(runif(50), nrow = 10, ncol = 5),
    frequencies = 1:5,
    times = seq(0, 0.9, by = 0.1)
  )
  p <- plotWaveletCoherence(result, title = "Custom Title",
                            fill_label = "Custom Label")
  expect_s3_class(p, "ggplot")
})

# ---- .convolve_mirror --------------------------------------------------------

test_that(".convolve_mirror preserves signal length", {
  x <- rnorm(100)
  kernel <- c(0.25, 0.5, 0.25)
  result <- PhysioCrossModal:::.convolve_mirror(x, kernel)
  expect_length(result, 100)
})

test_that(".convolve_mirror handles complex input", {
  x <- complex(real = rnorm(50), imaginary = rnorm(50))
  kernel <- c(0.25, 0.5, 0.25)
  result <- PhysioCrossModal:::.convolve_mirror(x, kernel)
  expect_length(result, 50)
  expect_true(is.complex(result))
})

test_that(".convolve_mirror with single-element kernel returns input", {
  x <- rnorm(50)
  result <- PhysioCrossModal:::.convolve_mirror(x, 1)
  expect_equal(result, x)
})
