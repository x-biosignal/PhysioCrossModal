# Multi-block latent fusion (JIVE, MFA), verified on synthetic blocks with a
# known shared + block-specific structure.

make_blocks <- function(n = 60, seed = 1) {
  set.seed(seed)
  shared <- matrix(rnorm(n * 2), n, 2)                 # joint structure (both blocks)
  ind1   <- matrix(rnorm(n * 1), n, 1)                 # kin-specific
  ind2   <- matrix(rnorm(n * 1), n, 1)                 # emg-specific
  kin <- shared %*% matrix(rnorm(2 * 8), 2, 8) +
         ind1  %*% matrix(rnorm(1 * 8), 1, 8) + matrix(rnorm(n * 8, 0, 0.15), n, 8)
  emg <- shared %*% matrix(rnorm(2 * 6), 2, 6) +
         ind2  %*% matrix(rnorm(1 * 6), 1, 6) + matrix(rnorm(n * 6, 0, 0.15), n, 6)
  list(blocks = list(kin = kin, emg = emg), shared = shared, ind1 = ind1, ind2 = ind2)
}

test_that("JIVE separates joint from individual variation", {
  d <- make_blocks()
  fit <- jive(d$blocks, rank_joint = 2, rank_individual = 1)
  expect_s3_class(fit, "jive")
  # both blocks carry substantial joint variance
  expect_true(all(fit$variance$joint > 0.3))
  # joint + individual + residual = 1 per block
  expect_equal(fit$variance$joint + fit$variance$individual + fit$variance$residual,
               rep(1, 2), tolerance = 1e-8)
  # joint scores recover the true shared subject scores (canonical correlation ~ 1)
  cc <- cancor(fit$joint_scores, d$shared)$cor
  expect_gt(min(cc), 0.95)
})

test_that("JIVE joint variance vanishes when blocks share nothing", {
  set.seed(2); n <- 60
  b1 <- matrix(rnorm(n * 6), n, 6); b2 <- matrix(rnorm(n * 6), n, 6)   # independent
  fit <- jive(list(a = b1, b = b2), rank_joint = 1, rank_individual = 1)
  # only spurious finite-sample joint variance (well below the >0.3 of the
  # genuinely-shared case)
  expect_lt(max(fit$variance$joint), 0.25)
})

test_that("JIVE individual structure maps to the right block", {
  # block 1 has an extra strong idiosyncratic component, block 2 does not
  set.seed(3); n <- 80
  shared <- matrix(rnorm(n * 1), n, 1)
  b1 <- shared %*% matrix(rnorm(5), 1, 5) +
        matrix(rnorm(n), n, 1) %*% matrix(rnorm(5), 1, 5) * 3 +   # big individual
        matrix(rnorm(n * 5, 0, 0.1), n, 5)
  b2 <- shared %*% matrix(rnorm(5), 1, 5) + matrix(rnorm(n * 5, 0, 0.1), n, 5)
  fit <- jive(list(b1 = b1, b2 = b2), rank_joint = 1, rank_individual = 1)
  expect_gt(fit$variance$individual[1], fit$variance$individual[2])
  expect_output(print(fit), "JIVE")
})

test_that("MFA global factor recovers the shared structure and both blocks load", {
  d <- make_blocks(seed = 4)
  m <- multipleFactorAnalysis(d$blocks, ncomp = 3)
  expect_s3_class(m, "mfa")
  expect_length(m$explained, 3)
  # first global factor scores relate to the shared structure
  cc <- cancor(m$scores[, 1:2], d$shared)$cor
  expect_gt(max(cc), 0.9)
  # both blocks contribute to the leading factor (neither dominates)
  expect_true(all(m$contributions[, 1] > 0.1))
})

test_that("MFA block weighting stops a high-variance block dominating", {
  set.seed(5); n <- 50
  s <- matrix(rnorm(n * 1), n, 1)
  small <- s %*% matrix(rnorm(4), 1, 4) + matrix(rnorm(n * 4, 0, 0.2), n, 4)
  big   <- (s %*% matrix(rnorm(4), 1, 4) + matrix(rnorm(n * 4, 0, 0.2), n, 4)) * 100
  m <- multipleFactorAnalysis(list(small = small, big = big), ncomp = 2)
  # despite the 100x scale, both blocks contribute comparably to factor 1
  expect_gt(min(m$contributions[, 1]), 0.2)
})

test_that("inputs must share the observation dimension", {
  expect_error(jive(list(matrix(0, 10, 3), matrix(0, 8, 3)), rank_joint = 1),
               "same number of rows")
  expect_error(multipleFactorAnalysis(list(matrix(0, 10, 3))), ">= 2")
})
