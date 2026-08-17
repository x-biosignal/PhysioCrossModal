# Coupled NMF, verified on non-negative blocks sharing a common drive.

test_that("coupled NMF recovers a shared latent drive from both modalities", {
  set.seed(1); n <- 80
  W <- matrix(runif(n * 2), n, 2)                      # shared temporal drive
  X <- W %*% matrix(runif(2 * 6), 2, 6) + matrix(runif(n * 6, 0, 0.02), n, 6)
  Y <- W %*% matrix(runif(2 * 4), 2, 4) + matrix(runif(n * 4, 0, 0.02), n, 4)
  fit <- coupledNMF(X, Y, rank = 2, n_restart = 4, seed = 1)
  expect_s3_class(fit, "coupled_nmf")
  expect_gt(fit$vaf_x, 0.95); expect_gt(fit$vaf_y, 0.95)
  # the recovered shared factor spans the true shared drive
  expect_gt(min(cancor(fit$shared, W)$cor), 0.9)
  expect_true(all(fit$shared >= 0) && all(fit$Hx >= 0) && all(fit$Hy >= 0))
})

test_that("combined VAF is high and dimensions are correct", {
  set.seed(2); n <- 60
  W <- matrix(runif(n * 3), n, 3)
  X <- W %*% matrix(runif(3 * 8), 3, 8); Y <- W %*% matrix(runif(3 * 5), 3, 5)
  fit <- coupledNMF(X, Y, rank = 3, n_restart = 3, seed = 2)
  expect_gt(fit$vaf, 0.98)
  expect_equal(dim(fit$shared), c(n, 3))
  expect_equal(dim(fit$Hx), c(3, 8)); expect_equal(dim(fit$Hy), c(3, 5))
  expect_output(print(fit), "Coupled NMF")
})

test_that("higher rank explains more variance", {
  set.seed(3); n <- 60
  W <- matrix(runif(n * 3), n, 3)
  X <- W %*% matrix(runif(3 * 6), 3, 6) + matrix(runif(n * 6, 0, 0.05), n, 6)
  Y <- W %*% matrix(runif(3 * 5), 3, 5) + matrix(runif(n * 5, 0, 0.05), n, 5)
  v1 <- coupledNMF(X, Y, rank = 1, n_restart = 3, seed = 3)$vaf
  v3 <- coupledNMF(X, Y, rank = 3, n_restart = 3, seed = 3)$vaf
  expect_gt(v3, v1)
})

test_that("coupled NMF rejects negative input and mismatched rows", {
  expect_error(coupledNMF(matrix(rnorm(30), 10, 3), matrix(runif(30), 10, 3), rank = 1),
               "non-negative")
  expect_error(coupledNMF(matrix(runif(30), 10, 3), matrix(runif(24), 8, 3), rank = 1),
               "same number of rows")
})
