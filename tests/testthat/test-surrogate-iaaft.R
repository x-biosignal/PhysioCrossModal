test_that("IAAFT surrogate preserves the power spectrum (< 1% MAE)", {
  set.seed(1)
  n <- 1024; t <- (0:(n - 1)) / 100
  x <- sin(2 * pi * 5 * t) + 0.5 * sin(2 * pi * 11 * t) + rnorm(n, 0, 0.3)
  target_amp <- Mod(stats::fft(x))
  maes <- vapply(seq_len(20), function(k) {
    s <- PhysioCrossModal:::.iaaft_surrogate(x)
    mean(abs(Mod(stats::fft(s)) - target_amp)) / mean(target_amp)
  }, numeric(1))
  expect_lt(max(maes), 0.01)
})

test_that("IAAFT surrogate preserves the value distribution (KS p > 0.1)", {
  set.seed(2)
  x <- rgamma(1000, shape = 2)          # non-Gaussian distribution
  s <- PhysioCrossModal:::.iaaft_surrogate(x)
  # the surrogate is a permutation of x -> exactly the same distribution
  expect_equal(sort(s), sort(x))
  expect_gt(suppressWarnings(stats::ks.test(s, x)$p.value), 0.1)
})

test_that("AAFT surrogate returns a permutation of the input", {
  set.seed(3)
  x <- rnorm(512) + sin(2 * pi * 3 * (0:511) / 100)
  s <- PhysioCrossModal:::.aaft_surrogate(x)
  expect_length(s, length(x))
  expect_equal(sort(s), sort(x))
})

test_that("surrogateTest with IAAFT: coupled significant, independent not", {
  sr <- 100; n <- 1024; t <- (0:(n - 1)) / sr

  set.seed(7)
  xc <- sin(2 * pi * 8 * t)
  yc <- sin(2 * pi * 8 * t + 0.7) + rnorm(n, 0, 0.4)   # genuinely coupled
  rc <- surrogateTest(xc, yc, sr = sr, method = "plv", freq_band = c(6, 10),
                      surrogate_type = "iaaft", n_surrogates = 200)
  expect_lt(rc$p_value, 0.05)

  set.seed(11)
  xi <- rnorm(n); yi <- rnorm(n)                         # independent
  ri <- surrogateTest(xi, yi, sr = sr, method = "plv", freq_band = c(6, 10),
                      surrogate_type = "iaaft", n_surrogates = 200)
  expect_gt(ri$p_value, 0.05)
})

test_that("surrogateTest and surrogateMatrixTest accept iaaft and aaft", {
  sr <- 100; n <- 512; t <- (0:(n - 1)) / sr
  set.seed(5)
  x <- sin(2 * pi * 8 * t); y <- sin(2 * pi * 8 * t + 0.5) + rnorm(n, 0, 0.4)

  for (st in c("iaaft", "aaft")) {
    res <- surrogateTest(x, y, sr = sr, method = "plv", freq_band = c(6, 10),
                         surrogate_type = st, n_surrogates = 50)
    expect_true(is.numeric(res$p_value))
    expect_true(res$p_value >= 0 && res$p_value <= 1)
  }

  pex <- PhysioExperiment(
    assays = list(raw = cbind(x, rnorm(n))),
    colData = S4Vectors::DataFrame(label = c("A", "B"), type = c("EEG", "EEG")),
    samplingRate = sr)
  pey <- PhysioExperiment(
    assays = list(raw = cbind(y, rnorm(n))),
    colData = S4Vectors::DataFrame(label = c("C", "D"), type = c("EEG", "EEG")),
    samplingRate = sr)
  mt <- surrogateMatrixTest(pex, pey, method = "plv", freq_band = c(6, 10),
                            surrogate_type = "iaaft", n_surrogates = 30)
  expect_type(mt, "list")
})
