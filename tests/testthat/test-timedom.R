# Tests for time-domain coupling functions (coupling-timedom.R)

# ---- crossCorrelation -------------------------------------------------------

test_that("crossCorrelation detects known lag in directed signals", {
  # y is a delayed copy of x with lag=5 at sr=500
  # make_directed_signals: y[i] = coupling * x[i - lag_samples]
  # In ccf(x, y) convention, lag h computes cor(x[t+h], y[t]).
  # Maximum occurs at h = -lag_samples (x leads y by lag_samples).
  sigs <- make_directed_signals(n = 5000, sr = 500, lag_samples = 5,
                                 coupling = 1.0, seed = 42)
  result <- crossCorrelation(sigs$x, sigs$y, sr = sigs$sr)

  expect_type(result, "list")
  expect_true("correlation" %in% names(result))
  expect_true("lags" %in% names(result))
  expect_true("lag_seconds" %in% names(result))
  expect_true("peak_lag" %in% names(result))
  expect_true("peak_lag_seconds" %in% names(result))
  expect_true("peak_correlation" %in% names(result))

  # Peak lag should be -5 samples (negative = x leads y)
  expect_equal(result$peak_lag, -5L)

  # Peak lag in seconds: -5 / 500 = -0.01
  expect_equal(result$peak_lag_seconds, -0.01, tolerance = 1e-6)
})

test_that("crossCorrelation of identical signals gives peak_lag=0 and r~1", {
  set.seed(123)
  x <- rnorm(2000)
  result <- crossCorrelation(x, x, sr = 500)

  expect_equal(result$peak_lag, 0L)
  expect_equal(result$peak_correlation, 1.0, tolerance = 1e-10)
})

test_that("crossCorrelation of independent signals has low peak correlation", {
  set.seed(100)
  x <- rnorm(5000)
  y <- rnorm(5000)
  result <- crossCorrelation(x, y, sr = 500)

  expect_true(abs(result$peak_correlation) < 0.3,
    info = paste("Peak correlation:", result$peak_correlation,
                 "(should be < 0.3 for independent signals)"))
})

test_that("max_lag parameter limits the lag range", {
  set.seed(42)
  x <- rnorm(1000)
  y <- rnorm(1000)

  result_default <- crossCorrelation(x, y, sr = 500)
  result_limited <- crossCorrelation(x, y, sr = 500, max_lag = 50L)

  # Default max_lag should be floor(1000/4) = 250
  expect_equal(max(result_default$lags), 250L)
  expect_equal(min(result_default$lags), -250L)

  # Limited max_lag should be 50
  expect_equal(max(result_limited$lags), 50L)
  expect_equal(min(result_limited$lags), -50L)

  # The limited result should have fewer lag values
  expect_true(length(result_limited$lags) < length(result_default$lags))
  expect_equal(length(result_limited$lags), 101L)  # -50 to 50 inclusive
})

test_that("normalize=FALSE gives unnormalized values", {
  set.seed(42)
  # Use signals with known variance for comparison
  x <- rnorm(2000) * 5  # larger variance
  y <- rnorm(2000) * 5

  result_norm <- crossCorrelation(x, y, sr = 500, normalize = TRUE)
  result_unnorm <- crossCorrelation(x, y, sr = 500, normalize = FALSE)

  # Normalized values should be bounded by [-1, 1] (approximately)
  expect_true(all(abs(result_norm$correlation) <= 1 + 1e-10))

  # Unnormalized values (covariance) can exceed [-1, 1] if variance > 1
  # The two results should have the same lag structure

  expect_equal(result_norm$lags, result_unnorm$lags)

  # Unnormalized and normalized should NOT be identical for non-unit-variance
  expect_false(isTRUE(all.equal(result_norm$correlation,
                                 result_unnorm$correlation)))
})

test_that("crossCorrelation works with numeric vectors directly", {
  set.seed(42)
  x <- rnorm(500)
  # y is a delayed copy of x: y[t] ~ x[t-3], so x leads y by 3 samples
  y <- c(rep(0, 3), x[1:497])

  result <- crossCorrelation(x, y, sr = 100)

  expect_type(result, "list")
  # ccf convention: peak at negative lag when x leads y
  expect_equal(result$peak_lag, -3L)
  expect_equal(result$peak_lag_seconds, -0.03, tolerance = 1e-6)
})

test_that("crossCorrelation returns consistent vector lengths", {
  set.seed(42)
  x <- rnorm(1000)
  y <- rnorm(1000)
  max_lag <- 100L

  result <- crossCorrelation(x, y, sr = 500, max_lag = max_lag)

  # All vectors should have the same length: 2 * max_lag + 1
  n_expected <- 2L * max_lag + 1L
  expect_length(result$correlation, n_expected)
  expect_length(result$lags, n_expected)
  expect_length(result$lag_seconds, n_expected)
})

test_that("crossCorrelation requires sr for numeric vectors", {
  expect_error(crossCorrelation(rnorm(100), rnorm(100)),
               "sr.*required")
})

# ---- slidingCrossCorrelation ------------------------------------------------

test_that("slidingCrossCorrelation produces correct number of windows", {
  set.seed(42)
  sr <- 500
  n_sec <- 10
  n <- sr * n_sec  # 5000 samples
  x <- rnorm(n)
  y <- rnorm(n)

  window_sec <- 1
  step_sec <- 0.5

  result <- slidingCrossCorrelation(x, y, sr = sr,
                                     window_sec = window_sec,
                                     step_sec = step_sec)

  # Expected number of windows: floor((n - window_samples) / step_samples) + 1
  window_samples <- window_sec * sr  # 500
  step_samples <- step_sec * sr      # 250
  expected_windows <- floor((n - window_samples) / step_samples) + 1L

  expect_type(result, "list")
  expect_true("correlations" %in% names(result))
  expect_true("times" %in% names(result))
  expect_true("lags" %in% names(result))
  expect_true("peak_lags" %in% names(result))
  expect_true("peak_correlations" %in% names(result))

  expect_equal(nrow(result$correlations), expected_windows)
  expect_length(result$times, expected_windows)
  expect_length(result$peak_lags, expected_windows)
  expect_length(result$peak_correlations, expected_windows)

  # Correlations matrix should have n_lags columns
  max_lag_expected <- floor(window_samples / 4L)
  n_lags_expected <- 2L * max_lag_expected + 1L
  expect_equal(ncol(result$correlations), n_lags_expected)
  expect_length(result$lags, n_lags_expected)
})

test_that("slidingCrossCorrelation peak_lags are consistent for stationary coupling", {
  # Create a stationary coupling: y is always a delayed copy of x
  # y[t] = x[t - lag_samples], so x leads y
  sr <- 500
  n <- 5000
  lag_samples <- 10L

  set.seed(42)
  x <- rnorm(n)
  y <- c(rep(0, lag_samples), x[1:(n - lag_samples)])

  result <- slidingCrossCorrelation(x, y, sr = sr,
                                     window_sec = 1, step_sec = 0.5)

  # Most windows (excluding edge effects) should detect the correct lag
  # ccf convention: peak at -lag_samples when x leads y
  n_windows <- length(result$peak_lags)
  if (n_windows > 2) {
    # Check inner windows (skip first and last which may have edge effects)
    inner_lags <- result$peak_lags[2:(n_windows - 1)]
    fraction_correct <- mean(inner_lags == -lag_samples)
    expect_true(fraction_correct > 0.8,
      info = paste("Fraction of windows with correct lag:",
                   round(fraction_correct, 3)))
  }
})

test_that("slidingCrossCorrelation with custom max_lag", {
  set.seed(42)
  sr <- 500
  x <- rnorm(5000)
  y <- rnorm(5000)

  result <- slidingCrossCorrelation(x, y, sr = sr,
                                     window_sec = 1, step_sec = 0.5,
                                     max_lag = 20L)

  # Lags should range from -20 to 20
  expect_equal(min(result$lags), -20L)
  expect_equal(max(result$lags), 20L)
  expect_equal(ncol(result$correlations), 41L)
})

test_that("slidingCrossCorrelation times are window centres", {
  sr <- 500
  n <- 5000
  set.seed(42)
  x <- rnorm(n)
  y <- rnorm(n)

  result <- slidingCrossCorrelation(x, y, sr = sr,
                                     window_sec = 2, step_sec = 1)

  # First window: samples 1..1000, centre at sample 500 => 499.5/500 ~ 0.999
  # The centre time should be approximately (window_sec / 2) for the first window
  # and increase by step_sec for each subsequent window
  expect_true(length(result$times) > 1)

  # Check that times are monotonically increasing
  diffs <- diff(result$times)
  expect_true(all(diffs > 0),
    info = "Window centre times should be monotonically increasing")

  # Steps between consecutive centres should be approximately step_sec
  expect_equal(diffs, rep(1, length(diffs)), tolerance = 0.01)
})

test_that("slidingCrossCorrelation rejects window larger than signal", {
  expect_error(
    slidingCrossCorrelation(rnorm(100), rnorm(100), sr = 500,
                             window_sec = 5),
    "exceeds"
  )
})

test_that("slidingCrossCorrelation works with numeric vectors", {
  set.seed(42)
  x <- rnorm(2000)
  y <- rnorm(2000)

  result <- slidingCrossCorrelation(x, y, sr = 500,
                                     window_sec = 1, step_sec = 0.5)

  expect_type(result, "list")
  expect_true(is.matrix(result$correlations))
  expect_true(nrow(result$correlations) > 0)
})
