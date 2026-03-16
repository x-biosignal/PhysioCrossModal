# Tests for coupling-wrapper.R and vis-coupling.R

# ---- couplingAnalysis dispatch -----------------------------------------------

test_that("couplingAnalysis dispatches to coherence", {
  sigs <- make_coupled_signals(
    sr1 = 500, sr2 = 500, n_sec = 10,
    coupling_freq = 20, coupling_strength = 0.8, noise_level = 0.2
  )

  result <- couplingAnalysis(sigs$x, sigs$y, method = "coherence",
                             sr = sigs$sr_x, nperseg = 256L)

  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true("confidence_limit" %in% names(result))
  expect_true("n_segments" %in% names(result))

  # Check that the coupling frequency shows high coherence

  idx_20 <- which.min(abs(result$frequencies - 20))
  expect_true(result$coherence[idx_20] > 0.5)
})

test_that("couplingAnalysis dispatches to plv", {
  sr <- 500
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t)
  y <- sin(2 * pi * 10 * t + pi / 4)

  result <- couplingAnalysis(x, y, method = "plv", sr = sr,
                             freq_band = c(8, 12))

  expect_type(result, "list")
  expect_true("plv" %in% names(result))
  expect_true("phase_diff" %in% names(result))
  expect_true(result$plv > 0.8)
})

test_that("couplingAnalysis dispatches to pli", {
  sr <- 500
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t)
  y <- sin(2 * pi * 10 * t + pi / 3)

  result <- couplingAnalysis(x, y, method = "pli", sr = sr,
                             freq_band = c(8, 12))

  expect_type(result, "list")
  expect_true("pli" %in% names(result))
  expect_true("phase_diff" %in% names(result))
})

test_that("couplingAnalysis dispatches to wpli", {
  sr <- 500
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t)
  y <- sin(2 * pi * 10 * t + pi / 4)

  result <- couplingAnalysis(x, y, method = "wpli", sr = sr,
                             freq_band = c(8, 12))

  expect_type(result, "list")
  expect_true("wpli" %in% names(result))
  expect_true("wpli_debiased" %in% names(result))
  expect_true("n_samples" %in% names(result))
})

test_that("couplingAnalysis dispatches to granger", {
  dsig <- make_directed_signals(n = 5000, sr = 500, lag_samples = 5,
                                 coupling = 0.6, seed = 42)

  result <- couplingAnalysis(dsig$x, dsig$y, method = "granger",
                             sr = dsig$sr, order = 10L)

  expect_type(result, "list")
  expect_true("gc_xy" %in% names(result))
  expect_true("gc_yx" %in% names(result))
  expect_true("net_gc" %in% names(result))
  expect_true("order" %in% names(result))
  # x drives y, so gc_xy should be positive
  expect_true(result$gc_xy > 0)
})

test_that("couplingAnalysis dispatches to crosscorrelation", {
  sr <- 500
  set.seed(1)
  x <- rnorm(1000)
  y <- c(rep(0, 10), x[1:990])

  result <- couplingAnalysis(x, y, method = "crosscorrelation", sr = sr)

  expect_type(result, "list")
  expect_true("correlation" %in% names(result))
  expect_true("lags" %in% names(result))
  expect_true("peak_lag" %in% names(result))
  expect_true("peak_correlation" %in% names(result))
})

test_that("couplingAnalysis dispatches to wavelet_coherence", {
  set.seed(42)
  sr <- 200
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t) + 0.2 * rnorm(n)

  result <- couplingAnalysis(x, y, method = "wavelet_coherence", sr = sr,
                             frequencies = seq(5, 15))

  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true("times" %in% names(result))
})

test_that("couplingAnalysis dispatches to wavelet_plv", {
  set.seed(42)
  sr <- 200
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- sin(2 * pi * 10 * t + pi / 4) + 0.2 * rnorm(n)

  result <- couplingAnalysis(x, y, method = "wavelet_plv", sr = sr,
                             frequencies = seq(5, 15))

  expect_type(result, "list")
  expect_true("plv" %in% names(result))
  expect_true("frequencies" %in% names(result))
})

test_that("couplingAnalysis errors on invalid method", {
  expect_error(
    couplingAnalysis(rnorm(100), rnorm(100), method = "invalid_method",
                     sr = 500)
  )
})

test_that("couplingAnalysis works with mpe argument", {
  pairs <- make_eeg_emg(n_sec = 10, sr_eeg = 500, sr_emg = 500)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- couplingAnalysis(
    mpe = mpe,
    modality_x = "EEG", modality_y = "EMG",
    channels_x = 1L, channels_y = 1L,
    method = "coherence",
    nperseg = 256L
  )

  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
  expect_true("frequencies" %in% names(result))
})


# ---- plotCouplingMatrix ------------------------------------------------------

test_that("plotCouplingMatrix returns a ggplot", {
  skip_if_not_installed("ggplot2")

  mat <- matrix(runif(16), nrow = 4,
                dimnames = list(paste0("Ch", 1:4), paste0("Ch", 1:4)))
  p <- plotCouplingMatrix(mat)
  expect_s3_class(p, "ggplot")
})

test_that("plotCouplingMatrix works without row/col names", {
  skip_if_not_installed("ggplot2")

  mat <- matrix(runif(9), nrow = 3)
  p <- plotCouplingMatrix(mat, title = "Test Matrix")
  expect_s3_class(p, "ggplot")
})

test_that("plotCouplingMatrix errors on non-matrix input", {
  skip_if_not_installed("ggplot2")

  expect_error(plotCouplingMatrix(1:10), "'mat' must be a numeric matrix")
})

# ---- plotCoherenceSpectrum ---------------------------------------------------

test_that("plotCoherenceSpectrum returns a ggplot", {
  skip_if_not_installed("ggplot2")

  sigs <- make_coupled_signals(sr1 = 500, sr2 = 500, n_sec = 10,
                                coupling_freq = 20, coupling_strength = 0.8,
                                noise_level = 0.2)
  coh_result <- coherence(sigs$x, sigs$y, sr = sigs$sr_x, nperseg = 256L)

  p <- plotCoherenceSpectrum(coh_result)
  expect_s3_class(p, "ggplot")
})

test_that("plotCoherenceSpectrum works without confidence limit", {
  skip_if_not_installed("ggplot2")

  fake_result <- list(
    coherence = runif(50),
    frequencies = seq(0, 250, length.out = 50)
  )
  p <- plotCoherenceSpectrum(fake_result, show_threshold = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plotCoherenceSpectrum errors on bad input", {
  skip_if_not_installed("ggplot2")

  expect_error(plotCoherenceSpectrum(list(a = 1)),
               "'result' must be a list with 'coherence' and 'frequencies'")
})

# ---- plotCouplingTimecourse --------------------------------------------------

test_that("plotCouplingTimecourse returns a ggplot", {
  skip_if_not_installed("ggplot2")

  sr <- 500
  set.seed(1)
  x <- rnorm(5000)
  y <- c(rep(0, 10), x[1:4990])
  slide_result <- slidingCrossCorrelation(x, y, sr = sr,
                                           window_sec = 1, step_sec = 0.5)

  p <- plotCouplingTimecourse(slide_result)
  expect_s3_class(p, "ggplot")
})

test_that("plotCouplingTimecourse errors on bad input", {
  skip_if_not_installed("ggplot2")

  expect_error(plotCouplingTimecourse(list(a = 1)),
               "'result' must be a list with 'times' and 'peak_correlations'")
})
