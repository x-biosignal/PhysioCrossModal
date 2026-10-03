# Two-block cross-modal methods (CCA, PLS), verified on blocks sharing a latent.

shared_blocks <- function(n = 60, k = 2, sd = 0.3, seed = 1) {
  set.seed(seed)
  z <- matrix(rnorm(n * k), n, k)                      # shared latent
  X <- z %*% matrix(rnorm(k * 6), k, 6) + matrix(rnorm(n * 6, 0, sd), n, 6)
  Y <- z %*% matrix(rnorm(k * 5), k, 5) + matrix(rnorm(n * 5, 0, sd), n, 5)
  list(X = X, Y = Y, z = z)
}

test_that("CCA recovers a strong cross-modal correlation for shared structure", {
  d <- shared_blocks()
  cc <- cca(d$X, d$Y)
  expect_s3_class(cc, "cca")
  expect_gt(cc$cor[1], 0.9)                            # top canonical correlation high
  expect_true(all(diff(cc$cor) <= 1e-8))              # decreasing
  # the top-2 canonical variates SPAN the 2-D shared latent subspace
  expect_gt(min(cancor(cc$xscores[, 1:2], d$z)$cor), 0.85)
})

test_that("CCA correlations are low for independent modalities", {
  set.seed(2); n <- 80
  X <- matrix(rnorm(n * 4), n, 4); Y <- matrix(rnorm(n * 4), n, 4)
  cc <- cca(X, Y)
  expect_lt(cc$cor[1], 0.6)                            # no genuine shared direction
})

test_that("CCA stays well-posed when features exceed observations (auto-ridge)", {
  set.seed(3); n <- 20
  z <- matrix(rnorm(n * 1), n, 1)
  X <- z %*% matrix(rnorm(30), 1, 30) + matrix(rnorm(n * 30, 0, .3), n, 30)  # p=30 > n
  Y <- z %*% matrix(rnorm(25), 1, 25) + matrix(rnorm(n * 25, 0, .3), n, 25)
  cc <- cca(X, Y)                                      # must not error / NaN
  expect_true(all(is.finite(cc$cor)))
  expect_true(all(cc$cor <= 1 + 1e-8))
  expect_gt(cc$cor[1], 0.5)
})

test_that("PLS finds the cross-modal covarying directions", {
  d <- shared_blocks(k = 1, seed = 4)
  pl <- plsBlocks(d$X, d$Y, ncomp = 2)
  expect_s3_class(pl, "pls_blocks")
  expect_gt(pl$correlation[1], 0.85)                  # first components strongly related
  expect_gt(pl$covariance[1], pl$covariance[2])       # covariance decreasing
  # the first PLS scores track the shared latent
  expect_gt(abs(cor(pl$xscores[, 1], d$z)), 0.85)
  expect_output(print(pl), "Two-block PLS")
})

test_that("PLS handles more features than observations", {
  set.seed(5); n <- 25
  z <- matrix(rnorm(n), n, 1)
  X <- z %*% matrix(rnorm(40), 1, 40) + matrix(rnorm(n * 40, 0, .3), n, 40)
  Y <- z %*% matrix(rnorm(35), 1, 35) + matrix(rnorm(n * 35, 0, .3), n, 35)
  pl <- plsBlocks(X, Y, ncomp = 1)
  expect_gt(pl$correlation[1], 0.8)
})

test_that("mismatched rows error", {
  expect_error(cca(matrix(0, 10, 3), matrix(0, 8, 3)), "same number of rows")
  expect_error(plsBlocks(matrix(0, 10, 3), matrix(0, 8, 3)), "same rows")
})
