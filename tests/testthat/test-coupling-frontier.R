band_noise <- function(n, seed = 1) {
  set.seed(seed)
  x <- as.numeric(stats::filter(rnorm(n), rep(1 / 5, 5), sides = 2))
  x[is.na(x)] <- 0
  x
}

# --- Phase Slope Index ----------------------------------------------------

test_that("PSI is directional, antisymmetric, and ~0 for independent signals", {
  sr <- 200; n <- 6000
  x <- band_noise(n)
  y <- c(rep(0, 8), x[seq_len(n - 8)])                # y lags x (x leads)
  z <- band_noise(n, seed = 99)                       # independent

  psi_xy <- phaseSlopeIndex(x, y, sr = sr, freq_band = c(1, 40))
  psi_yx <- phaseSlopeIndex(y, x, sr = sr, freq_band = c(1, 40))
  psi_ind <- phaseSlopeIndex(x, z, sr = sr, freq_band = c(1, 40))

  expect_gt(psi_xy, 0)                                # x leads y -> positive
  expect_equal(psi_xy, -psi_yx, tolerance = 1e-6)     # antisymmetric
  expect_lt(abs(psi_ind), abs(psi_xy) / 2)            # independent much smaller
  expect_true(abs(psi_xy) <= 1 + 1e-9)                # normalized range
})

# --- Orthogonalised AEC (leakage-corrected) -------------------------------

test_that("orthogonalized AEC removes pure zero-lag leakage", {
  sr <- 200; n <- 6000
  x <- band_noise(n)
  y_leak <- 0.8 * x                                   # pure instantaneous mixing
  aec_leak <- orthogonalizedAEC(x, y_leak, sr = sr, freq_band = c(8, 40))
  expect_lt(abs(aec_leak), 0.2)                       # leakage suppressed

  expect_true(aec_leak >= -1 && aec_leak <= 1)
})

test_that("orthogonalized AEC detects a genuine lagged envelope coupling", {
  sr <- 200; n <- 8000; t <- seq_len(n) / sr
  set.seed(3)
  env <- 1 + 0.8 * sin(2 * pi * 0.1 * t)              # shared slow envelope
  carrier <- function(ph) sin(2 * pi * 20 * t + ph)
  x <- env * carrier(0) + rnorm(n, sd = 0.2)
  y <- env * carrier(pi / 3) + rnorm(n, sd = 0.2)     # same envelope, phase-lagged
  aec <- orthogonalizedAEC(x, y, sr = sr, freq_band = c(15, 25))
  expect_gt(aec, 0.3)
})

# --- ciPLV ----------------------------------------------------------------

test_that("ciPLV is ~0 for zero-lag and high for a constant phase lag", {
  sr <- 200; n <- 6000; t <- seq_len(n) / sr
  set.seed(4)
  base <- sin(2 * pi * 10 * t)
  x <- base + rnorm(n, sd = 0.2)
  y0 <- base + rnorm(n, sd = 0.2)                     # zero-lag
  yl <- sin(2 * pi * 10 * t + pi / 2) + rnorm(n, sd = 0.2)  # 90-deg lag

  expect_lt(ciPLV(x, y0, sr = sr, freq_band = c(8, 12)), 0.3)
  expect_gt(ciPLV(x, yl, sr = sr, freq_band = c(8, 12)), 0.6)
})

# --- Pairwise phase consistency -------------------------------------------

test_that("PPC ~1 for locked, ~0 for independent phases", {
  sr <- 200; n <- 6000; t <- seq_len(n) / sr
  set.seed(5)
  x <- sin(2 * pi * 10 * t)
  y_lock <- sin(2 * pi * 10 * t + pi / 4) + rnorm(n, sd = 0.1)
  y_ind  <- sin(2 * pi * 6.3 * t) + rnorm(n, sd = 0.1)

  expect_gt(pairwisePhaseConsistency(x, y_lock, sr = sr, freq_band = c(8, 12)), 0.7)
  expect_lt(pairwisePhaseConsistency(x, y_ind, sr = sr, freq_band = c(8, 12)), 0.3)
})

# --- LEiDA dynamic states -------------------------------------------------

test_that("LEiDA recovers two connectivity regimes", {
  sr <- 100; n <- 3000; t <- seq_len(n) / sr; half <- n / 2
  set.seed(6)
  a <- sin(2 * pi * 10 * t)
  b <- sin(2 * pi * 10 * t + 0.2)
  # regime 1 (first half): ch1-2 coherent; regime 2: ch3-4 coherent
  X <- cbind(
    a + rnorm(n, sd = 0.3),
    c(a[1:half], rnorm(half)) + rnorm(n, sd = 0.3),
    c(rnorm(half), b[(half + 1):n]) + rnorm(n, sd = 0.3),
    b + rnorm(n, sd = 0.3))
  res <- leidaStates(X, sr = sr, freq_band = c(8, 12), n_states = 2, seed = 1)

  expect_s3_class(res, "leida")
  expect_equal(length(res$states), n)
  expect_equal(sum(res$occupancy), 1, tolerance = 1e-8)
  expect_equal(dim(res$centroids), c(2, 4))
  # the dominant state differs between the two halves
  first_mode  <- which.max(tabulate(res$states[1:half], 2))
  second_mode <- which.max(tabulate(res$states[(half + 1):n], 2))
  expect_false(first_mode == second_mode)
})
