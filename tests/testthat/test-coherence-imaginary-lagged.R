test_that("zero-lag identical signals: imaginary ~0, magnitude ~1", {
  sr <- 200; t <- seq(0, 20 - 1 / sr, by = 1 / sr); n <- length(t)
  set.seed(1)
  x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(n)
  y <- x                                                 # zero-lag identical

  mg <- coherence(x, y, sr = sr, freq_range = c(8, 12), type = "magnitude")
  im <- coherence(x, y, sr = sr, freq_range = c(8, 12), type = "imaginary")

  expect_gt(mean(mg$coherence), 0.99)                    # magnitude ~ 1
  expect_lt(max(abs(im$coherence)), 1e-3)                # imaginary ~ 0
})

test_that("90-degree phase-shifted narrowband pair: imaginary near max in band", {
  sr <- 200; t <- seq(0, 20 - 1 / sr, by = 1 / sr); n <- length(t)
  set.seed(2)
  a <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  b <- sin(2 * pi * 10 * t + pi / 2) + 0.2 * rnorm(n)    # 90-degree lag

  im <- coherence(a, b, sr = sr, freq_range = c(5, 15), type = "imaginary")
  peak_freq <- im$frequencies[which.max(abs(im$coherence))]
  expect_gt(max(abs(im$coherence)), 0.9)                 # near maximal
  expect_lt(abs(peak_freq - 10), 1)                      # in the 10 Hz band
})

test_that("lagged coherence matches the Pascual-Marqui closed form (< 1e-6)", {
  sr <- 200; t <- seq(0, 20 - 1 / sr, by = 1 / sr); n <- length(t)
  set.seed(3)
  a <- sin(2 * pi * 10 * t) + 0.3 * rnorm(n)
  b <- sin(2 * pi * 10 * t + 0.6) + 0.3 * rnorm(n)

  mg <- coherence(a, b, sr = sr, type = "magnitude")$coherence   # |C|^2
  im <- coherence(a, b, sr = sr, type = "imaginary")$coherence   # Im(C)
  lg <- coherence(a, b, sr = sr, type = "lagged")$coherence

  # |C|^2 = Re(C)^2 + Im(C)^2  =>  Re(C)^2 = |C|^2 - Im(C)^2
  re2 <- pmax(mg - im^2, 0)
  closed <- im / sqrt(pmax(1 - re2, .Machine$double.eps))
  expect_lt(max(abs(lg - closed)), 1e-6)
})

test_that("laggedCoherence wrapper returns lagged coherence", {
  sr <- 200; t <- seq(0, 10 - 1 / sr, by = 1 / sr); n <- length(t)
  x <- sin(2 * pi * 10 * t); y <- sin(2 * pi * 10 * t + pi / 2)
  lc <- laggedCoherence(x, y, sr = sr, freq_range = c(8, 12))
  expect_equal(lc$type, "lagged")
  expect_true(is.numeric(lc$coherence))
})

test_that("type flows through coherenceMatrix and couplingMatrix", {
  sr <- 200; t <- seq(0, 15 - 1 / sr, by = 1 / sr); n <- length(t)
  set.seed(4)
  a <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  b <- sin(2 * pi * 10 * t + pi / 2) + 0.2 * rnorm(n)
  pex <- PhysioExperiment(
    assays = list(raw = cbind(a, rnorm(n))),
    colData = S4Vectors::DataFrame(label = c("A", "B"), type = c("EEG", "EEG")),
    samplingRate = sr)
  pey <- PhysioExperiment(
    assays = list(raw = cbind(b, rnorm(n))),
    colData = S4Vectors::DataFrame(label = c("C", "D"), type = c("EEG", "EEG")),
    samplingRate = sr)

  cm <- coherenceMatrix(pex, pey, type = "imaginary")
  expect_equal(cm$spectra[[1]][[1]]$type, "imaginary")
  expect_true(cm$matrix[1, 1] >= 0)                      # abs peak magnitude

  cpm <- couplingMatrix(pex, pey, method = "coherence", type = "lagged")
  expect_true(is.matrix(cpm$matrix))
})

test_that("backward compatibility: default coherence is magnitude-squared", {
  sr <- 200; t <- seq(0, 10 - 1 / sr, by = 1 / sr); n <- length(t)
  set.seed(5)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  res <- coherence(x, y, sr = sr)
  expect_equal(res$type, "magnitude")
  expect_true(all(res$coherence >= 0 & res$coherence <= 1))
})

test_that("plotCoherenceSpectrum renders for imaginary/lagged types", {
  skip_if_not_installed("ggplot2")
  sr <- 200; t <- seq(0, 10 - 1 / sr, by = 1 / sr); n <- length(t)
  x <- sin(2 * pi * 10 * t); y <- sin(2 * pi * 10 * t + pi / 2)
  for (ty in c("imaginary", "lagged")) {
    res <- coherence(x, y, sr = sr, freq_range = c(5, 15), type = ty)
    p <- plotCoherenceSpectrum(res)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})
