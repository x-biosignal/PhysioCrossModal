# Tests for directed coupling: Granger causality (coupling-directed.R)

# ---- Basic directional coupling detection ------------------------------------

test_that("x drives y with lag=5: GC(x->y) > GC(y->x)", {
  sigs <- make_directed_signals(n = 5000, sr = 500, lag_samples = 5,
                                 coupling = 0.6, seed = 42)
  result <- grangerCausality(sigs$x, sigs$y, sr = sigs$sr, order = 10L)

  expect_type(result, "list")
  expect_true("gc_xy" %in% names(result))
  expect_true("gc_yx" %in% names(result))
  expect_true("net_gc" %in% names(result))
  expect_true("order" %in% names(result))

  # x drives y, so gc_xy should be substantially larger than gc_yx

  expect_true(result$gc_xy > result$gc_yx,
    info = sprintf("gc_xy = %.4f, gc_yx = %.4f", result$gc_xy, result$gc_yx))
})

test_that("x drives y: net_gc > 0", {
  sigs <- make_directed_signals(n = 5000, sr = 500, lag_samples = 5,
                                 coupling = 0.6, seed = 42)
  result <- grangerCausality(sigs$x, sigs$y, sr = sigs$sr, order = 10L)

  # net_gc = gc_xy - gc_yx should be positive when x drives y
  expect_true(result$net_gc > 0,
    info = sprintf("net_gc = %.4f, expected > 0", result$net_gc))
})

# ---- Independent signals ----------------------------------------------------

test_that("independent signals: both GC values near 0", {
  set.seed(123)
  n <- 5000
  x <- rnorm(n)
  y <- rnorm(n)

  result <- grangerCausality(x, y, sr = 500, order = 5L)

  # For independent signals, both GC values should be near 0
  expect_true(result$gc_xy < 0.1,
    info = sprintf("gc_xy = %.4f for independent signals, expected < 0.1",
                   result$gc_xy))
  expect_true(result$gc_yx < 0.1,
    info = sprintf("gc_yx = %.4f for independent signals, expected < 0.1",
                   result$gc_yx))
})

# ---- Symmetric coupling -----------------------------------------------------

test_that("symmetric coupling: net_gc approximately 0", {
  set.seed(99)
  n <- 5000
  # Create mutually coupled signals with equal strength
  x <- numeric(n)
  y <- numeric(n)
  x[1:3] <- rnorm(3)
  y[1:3] <- rnorm(3)
  for (i in 4:n) {
    x[i] <- 0.4 * y[i - 2] + 0.6 * rnorm(1)
    y[i] <- 0.4 * x[i - 2] + 0.6 * rnorm(1)
  }

  result <- grangerCausality(x, y, sr = 500, order = 5L)

  # net_gc should be approximately 0 for symmetric coupling
  expect_true(abs(result$net_gc) < 0.3,
    info = sprintf("net_gc = %.4f for symmetric coupling, expected |net_gc| < 0.3",
                   result$net_gc))
})

# ---- Order parameter effects -------------------------------------------------

test_that("order parameter affects results (higher order for longer lags)", {
  # Signal with a longer lag (10 samples) -- low order can't capture it well
  sigs <- make_directed_signals(n = 5000, sr = 500, lag_samples = 10,
                                 coupling = 0.6, seed = 42)

  result_low <- grangerCausality(sigs$x, sigs$y, sr = sigs$sr, order = 3L)
  result_high <- grangerCausality(sigs$x, sigs$y, sr = sigs$sr, order = 15L)

  # Both should return valid results
  expect_type(result_low, "list")
  expect_type(result_high, "list")
  expect_equal(result_low$order, 3L)
  expect_equal(result_high$order, 15L)

  # Higher order should capture the lag=10 coupling better
  # With order=3, the model cannot represent a 10-sample lag, so gc_xy should
  # be smaller than with order=15
  expect_true(result_high$gc_xy > result_low$gc_xy,
    info = sprintf("high order gc_xy = %.4f, low order gc_xy = %.4f",
                   result_high$gc_xy, result_low$gc_xy))
})

# ---- Numeric vector input ----------------------------------------------------

test_that("works with numeric vectors", {
  set.seed(42)
  n <- 3000
  x <- rnorm(n)
  y <- numeric(n)
  for (i in 6:n) y[i] <- 0.5 * x[i - 5] + 0.5 * rnorm(1)

  result <- grangerCausality(x, y, sr = 250, order = 8L)

  expect_type(result, "list")
  expect_true(is.numeric(result$gc_xy))
  expect_true(is.numeric(result$gc_yx))
  expect_true(is.numeric(result$net_gc))
  expect_equal(result$order, 8L)
  expect_true(result$gc_xy > 0)
})

# ---- Spectral Granger causality ---------------------------------------------

test_that("spectral GC with freq_range detects directional coupling", {
  sigs <- make_directed_signals(n = 5000, sr = 500, lag_samples = 5,
                                 coupling = 0.6, seed = 42)

  result <- grangerCausality(sigs$x, sigs$y, sr = sigs$sr, order = 10L,
                             freq_range = c(1, 100))

  expect_type(result, "list")
  # Spectral GC should also detect the directionality
  expect_true(result$gc_xy > result$gc_yx,
    info = sprintf("spectral gc_xy = %.4f, gc_yx = %.4f",
                   result$gc_xy, result$gc_yx))
  expect_true(result$net_gc > 0)
})

# ---- Input validation --------------------------------------------------------

test_that("grangerCausality validates order parameter", {
  x <- rnorm(100)
  y <- rnorm(100)

  expect_error(grangerCausality(x, y, sr = 100, order = 0L), "positive")
  expect_error(grangerCausality(x, y, sr = 100, order = -1L), "positive")
})

test_that("grangerCausality validates freq_range", {
  x <- rnorm(500)
  y <- rnorm(500)

  expect_error(grangerCausality(x, y, sr = 100, freq_range = c(20, 10)),
               "freq_range")
  expect_error(grangerCausality(x, y, sr = 100, freq_range = c(10)),
               "freq_range")
})

test_that("grangerCausality rejects too-short signals", {
  x <- rnorm(10)
  y <- rnorm(10)

  expect_error(grangerCausality(x, y, sr = 100, order = 5L), "too short")
})

# ---- Return structure --------------------------------------------------------

test_that("grangerCausality returns correct structure", {
  sigs <- make_directed_signals(n = 2000, sr = 500, lag_samples = 5,
                                 coupling = 0.5, seed = 1)
  result <- grangerCausality(sigs$x, sigs$y, sr = sigs$sr, order = 5L)

  expect_type(result, "list")
  expect_named(result, c("gc_xy", "gc_yx", "net_gc", "order"))
  expect_true(is.numeric(result$gc_xy) && length(result$gc_xy) == 1)
  expect_true(is.numeric(result$gc_yx) && length(result$gc_yx) == 1)
  expect_true(is.numeric(result$net_gc) && length(result$net_gc) == 1)
  expect_true(is.integer(result$order) && result$order == 5L)

  # GC values should be non-negative
  expect_true(result$gc_xy >= 0)
  expect_true(result$gc_yx >= 0)

  # net_gc should equal gc_xy - gc_yx
  expect_equal(result$net_gc, result$gc_xy - result$gc_yx, tolerance = 1e-12)
})

# ---- Internal helpers: direct tests -----------------------------------------

test_that(".gc_time_domain detects known causal signal", {
  set.seed(42)
  n <- 3000
  driver <- rnorm(n)
  target <- numeric(n)
  for (i in 4:n) target[i] <- 0.7 * driver[i - 3] + 0.3 * rnorm(1)

  gc_val <- PhysioCrossModal:::.gc_time_domain(driver - mean(driver),
                                                target - mean(target),
                                                order = 5L)
  expect_true(gc_val > 0.1,
    info = paste("GC for causal signal:", round(gc_val, 4)))
})

test_that(".gc_time_domain near 0 for independent signals", {
  set.seed(42)
  n <- 3000
  driver <- rnorm(n)
  target <- rnorm(n)

  gc_val <- PhysioCrossModal:::.gc_time_domain(driver, target, order = 5L)
  expect_true(gc_val < 0.1,
    info = paste("GC for independent signals:", round(gc_val, 4)))
})

test_that(".ols_residual_var returns correct variance for simple case", {
  # y = 2*x + noise
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  noise <- rnorm(n, sd = 0.5)
  y <- 2 * x + noise

  X <- matrix(x, ncol = 1)
  var_resid <- PhysioCrossModal:::.ols_residual_var(X, y)

  # Residual variance should be close to noise variance
  expect_true(abs(var_resid - 0.25) < 0.15,
    info = paste("Residual var:", round(var_resid, 3), "expected ~0.25"))
})

test_that(".ols_coefficients returns correct coefficients for simple case", {
  # y = 3*x1 - 1*x2
  set.seed(42)
  n <- 1000
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 3 * x1 - 1 * x2 + rnorm(n, sd = 0.1)

  X <- cbind(x1, x2)
  beta <- PhysioCrossModal:::.ols_coefficients(X, y)

  expect_length(beta, 2)
  expect_true(abs(beta[1] - 3) < 0.1,
    info = paste("Beta1:", round(beta[1], 3), "expected ~3"))
  expect_true(abs(beta[2] - (-1)) < 0.1,
    info = paste("Beta2:", round(beta[2], 3), "expected ~-1"))
})

test_that(".ols_residual_var with near-singular matrix uses ridge penalty", {
  # Create nearly collinear predictors
  set.seed(42)
  n <- 100
  x1 <- rnorm(n)
  x2 <- x1 + rnorm(n, sd = 1e-8)  # nearly identical to x1
  y <- rnorm(n)

  X <- cbind(x1, x2)
  # Should not error despite near-singularity
  var_resid <- PhysioCrossModal:::.ols_residual_var(X, y)
  expect_true(is.finite(var_resid))
  expect_true(var_resid > 0)
})
