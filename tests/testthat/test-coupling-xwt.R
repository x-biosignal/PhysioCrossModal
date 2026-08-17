test_that("XWT concentrates power at a shared oscillation and is significant", {
  sr <- 100; t <- seq(0, 20, by = 1 / sr)
  set.seed(1)
  x <- sin(2 * pi * 10 * t) + rnorm(length(t), sd = 0.3)
  y <- sin(2 * pi * 10 * t - pi / 2) + rnorm(length(t), sd = 0.3)  # y lags x by 90 deg
  freqs <- seq(4, 16, by = 1)
  set.seed(10)
  xw <- crossWaveletTransform(x, y, sr = sr, frequencies = freqs,
                              n_surrogates = 60)

  expect_s3_class(xw, "cross_wavelet")
  expect_equal(dim(xw$power), c(length(t), length(freqs)))

  # mean power is maximal at the 10 Hz row
  band_power <- colMeans(xw$power)
  expect_equal(freqs[which.max(band_power)], 10)

  # significance is common at 10 Hz, rare off-band
  i10 <- which(freqs == 10); i5 <- which(freqs == 5)
  expect_gt(mean(xw$significant[, i10]), 0.7)
  expect_lt(mean(xw$significant[, i5]), 0.3)
})

test_that("XWT phase recovers the imposed lead/lag", {
  sr <- 100; t <- seq(0, 20, by = 1 / sr)
  set.seed(2)
  x <- sin(2 * pi * 10 * t)
  y <- sin(2 * pi * 10 * t - pi / 2)                 # y lags x by pi/2
  freqs <- seq(6, 14, by = 1)
  xw <- crossWaveletTransform(x, y, sr = sr, frequencies = freqs,
                              significance = FALSE)
  i10 <- which(freqs == 10)
  # away from edges, relative phase ~ +pi/2 (x leads y)
  mid <- seq(sr, length(t) - sr)
  ph <- median(xw$phase[mid, i10])
  expect_equal(ph, pi / 2, tolerance = 0.3)
})

test_that("two independent noise signals hit the ~5% false-positive rate", {
  sr <- 100; t <- seq(0, 20, by = 1 / sr)
  set.seed(3)
  x <- rnorm(length(t))                              # independent noise
  z <- rnorm(length(t))                              # independent noise
  set.seed(11)
  xw <- crossWaveletTransform(x, z, sr = sr, frequencies = seq(4, 16, by = 1),
                              n_surrogates = 60, siglvl = 0.95)
  # with no true coupling, significance should sit near the nominal 5%
  expect_lt(mean(xw$significant), 0.12)
})
