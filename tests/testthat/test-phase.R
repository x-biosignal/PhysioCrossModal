# Tests for phase coupling functions (coupling-phase.R)
#
# Tests PLV, PLI, and wPLI with known synthetic signals.

# ---- PLV: Phase Locking Value ------------------------------------------------

test_that("PLV of two identical-phase sinusoids is approximately 1.0", {
  sr <- 500
  n_sec <- 4
  n <- sr * n_sec
  t <- seq(0, n_sec - 1 / sr, length.out = n)
  freq <- 10

  x <- sin(2 * pi * freq * t)
  y <- sin(2 * pi * freq * t)  # same frequency, same phase

  result <- phaseLockingValue(x, y, sr = sr, freq_band = c(8, 12))

  expect_type(result, "list")

  expect_true("plv" %in% names(result))
  expect_true("phase_diff" %in% names(result))
  expect_length(result$phase_diff, n)

  # PLV should be very high (> 0.9) for identical sinusoids
  expect_true(result$plv > 0.9,
    info = paste("PLV of identical sinusoids:", round(result$plv, 4)))
})

test_that("PLV of two independent noise signals is low", {
  set.seed(123)
  sr <- 500
  n <- 5000

  x <- rnorm(n)
  y <- rnorm(n)

  result <- phaseLockingValue(x, y, sr = sr, freq_band = c(8, 12))

  expect_true(result$plv < 0.3,
    info = paste("PLV of independent noise:", round(result$plv, 4)))
})

test_that("PLV of phase-shifted sinusoids is still high", {
  sr <- 500
  n_sec <- 4
  n <- sr * n_sec
  t <- seq(0, n_sec - 1 / sr, length.out = n)
  freq <- 10

  x <- sin(2 * pi * freq * t)
  y <- sin(2 * pi * freq * t + pi / 4)  # constant phase offset

  result <- phaseLockingValue(x, y, sr = sr, freq_band = c(8, 12))

  # Phase-locked signals with constant phase offset should still have high PLV
  expect_true(result$plv > 0.9,
    info = paste("PLV with pi/4 offset:", round(result$plv, 4)))
})

# ---- PLI: Phase Lag Index ----------------------------------------------------

test_that("PLI of phase-locked signals with non-zero lag is high", {
  sr <- 500
  n_sec <- 4
  n <- sr * n_sec
  t <- seq(0, n_sec - 1 / sr, length.out = n)
  freq <- 10

  x <- sin(2 * pi * freq * t)
  y <- sin(2 * pi * freq * t + pi / 3)  # non-zero phase lag

  result <- phaseLagIndex(x, y, sr = sr, freq_band = c(8, 12))

  expect_type(result, "list")
  expect_true("pli" %in% names(result))
  expect_true("phase_diff" %in% names(result))
  expect_length(result$phase_diff, n)

  # Non-zero lag should give high PLI
  expect_true(result$pli > 0.5,
    info = paste("PLI with pi/3 lag:", round(result$pli, 4)))
})

test_that("PLI of zero-lag coupling is approximately 0 (symmetric)", {
  sr <- 500
  n_sec <- 4
  n <- sr * n_sec
  t <- seq(0, n_sec - 1 / sr, length.out = n)
  freq <- 10

  x <- sin(2 * pi * freq * t)
  y <- sin(2 * pi * freq * t)  # zero phase lag

  result <- phaseLagIndex(x, y, sr = sr, freq_band = c(8, 12))

  # Zero-lag coupling: sin(phase_diff) is ~0 for all t, so PLI should be ~0
  expect_true(result$pli < 0.3,
    info = paste("PLI of zero-lag coupling:", round(result$pli, 4)))
})

test_that("PLI values are in [0, 1]", {
  set.seed(42)
  sr <- 500
  n <- 2000

  x <- rnorm(n)
  y <- rnorm(n)

  result <- phaseLagIndex(x, y, sr = sr, freq_band = c(8, 12))

  expect_true(result$pli >= 0 && result$pli <= 1,
    info = paste("PLI value:", round(result$pli, 4)))
})

# ---- wPLI: Weighted Phase Lag Index ------------------------------------------

test_that("wPLI returns value in [0, 1]", {
  sr <- 500
  n_sec <- 4
  n <- sr * n_sec
  t <- seq(0, n_sec - 1 / sr, length.out = n)
  freq <- 10

  x <- sin(2 * pi * freq * t)
  y <- sin(2 * pi * freq * t + pi / 4)

  result <- weightedPLI(x, y, sr = sr, freq_band = c(8, 12))

  expect_type(result, "list")
  expect_true("wpli" %in% names(result))
  expect_true("wpli_debiased" %in% names(result))
  expect_true("n_samples" %in% names(result))

  expect_true(result$wpli >= 0 && result$wpli <= 1,
    info = paste("wPLI value:", round(result$wpli, 4)))
  expect_equal(result$n_samples, n)
})

test_that("wPLI debiased version returns a valid value", {
  sr <- 500
  n_sec <- 4
  n <- sr * n_sec
  t <- seq(0, n_sec - 1 / sr, length.out = n)
  freq <- 10

  x <- sin(2 * pi * freq * t)
  y <- sin(2 * pi * freq * t + pi / 3)

  result_debiased <- weightedPLI(x, y, sr = sr, freq_band = c(8, 12),
                                 debiased = TRUE)
  result_no_debias <- weightedPLI(x, y, sr = sr, freq_band = c(8, 12),
                                  debiased = FALSE)

  # Debiased result should be present when debiased=TRUE

  expect_true(!is.null(result_debiased$wpli_debiased))
  expect_true(is.numeric(result_debiased$wpli_debiased))

  # Debiased result should be NULL when debiased=FALSE
  expect_null(result_no_debias$wpli_debiased)

  # For large N and high coupling, debiased should be close to the raw value
  expect_true(is.finite(result_debiased$wpli_debiased))
})

test_that("wPLI of phase-locked signals is high", {
  sr <- 500
  n_sec <- 4
  n <- sr * n_sec
  t <- seq(0, n_sec - 1 / sr, length.out = n)
  freq <- 10

  x <- sin(2 * pi * freq * t)
  y <- sin(2 * pi * freq * t + pi / 4)

  result <- weightedPLI(x, y, sr = sr, freq_band = c(8, 12))

  expect_true(result$wpli > 0.5,
    info = paste("wPLI of phase-locked signals:", round(result$wpli, 4)))
})

# ---- Input validation --------------------------------------------------------

test_that("freq_band is required (missing causes error)", {
  sr <- 500
  x <- rnorm(1000)
  y <- rnorm(1000)

  expect_error(phaseLockingValue(x, y, sr = sr),
               "freq_band")
  expect_error(phaseLagIndex(x, y, sr = sr),
               "freq_band")
  expect_error(weightedPLI(x, y, sr = sr),
               "freq_band")
})

test_that("freq_band must be numeric vector of length 2", {
  sr <- 500
  x <- rnorm(1000)
  y <- rnorm(1000)

  expect_error(phaseLockingValue(x, y, sr = sr, freq_band = 10),
               "length 2")
  expect_error(phaseLagIndex(x, y, sr = sr, freq_band = c(1, 2, 3)),
               "length 2")
  expect_error(weightedPLI(x, y, sr = sr, freq_band = "alpha"),
               "numeric")
})

# ---- Works with two numeric vectors and sr -----------------------------------

test_that("All functions work with two numeric vectors and sr parameter", {
  set.seed(42)
  sr <- 250
  n_sec <- 4
  n <- sr * n_sec
  t <- seq(0, n_sec - 1 / sr, length.out = n)
  freq <- 10

  x <- sin(2 * pi * freq * t) + 0.1 * rnorm(n)
  y <- sin(2 * pi * freq * t + pi / 6) + 0.1 * rnorm(n)

  # PLV
  plv_result <- phaseLockingValue(x, y, sr = sr, freq_band = c(8, 12))
  expect_true(is.numeric(plv_result$plv))
  expect_true(plv_result$plv >= 0 && plv_result$plv <= 1)
  expect_length(plv_result$phase_diff, n)

  # PLI
  pli_result <- phaseLagIndex(x, y, sr = sr, freq_band = c(8, 12))
  expect_true(is.numeric(pli_result$pli))
  expect_true(pli_result$pli >= 0 && pli_result$pli <= 1)
  expect_length(pli_result$phase_diff, n)

  # wPLI
  wpli_result <- weightedPLI(x, y, sr = sr, freq_band = c(8, 12))
  expect_true(is.numeric(wpli_result$wpli))
  expect_true(wpli_result$wpli >= 0 && wpli_result$wpli <= 1)
  expect_true(is.integer(wpli_result$n_samples))
})

# ---- Edge cases and consistency ----------------------------------------------

test_that("PLV phase_diff values are within [-pi, pi]", {
  sr <- 500
  n <- 2000
  t <- seq(0, (n - 1) / sr, length.out = n)

  x <- sin(2 * pi * 10 * t)
  y <- sin(2 * pi * 10 * t + pi / 2)

  result <- phaseLockingValue(x, y, sr = sr, freq_band = c(8, 12))

  expect_true(all(result$phase_diff >= -pi - 1e-10))
  expect_true(all(result$phase_diff <= pi + 1e-10))
})

test_that("wavelet method is not yet implemented", {
  sr <- 500
  x <- rnorm(1000)
  y <- rnorm(1000)

  expect_error(phaseLockingValue(x, y, sr = sr, freq_band = c(8, 12),
                                 method = "wavelet"),
               "not yet implemented")
})
