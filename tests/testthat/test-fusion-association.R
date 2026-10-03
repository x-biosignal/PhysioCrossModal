# Cross-modal whole-block association (RV, distance correlation, RSA).

shared_XY <- function(n = 50, sd = 0.3, seed = 1) {
  set.seed(seed); z <- matrix(rnorm(n * 2), n, 2)
  X <- z %*% matrix(rnorm(2 * 5), 2, 5) + matrix(rnorm(n * 5, 0, sd), n, 5)
  Y <- z %*% matrix(rnorm(2 * 4), 2, 4) + matrix(rnorm(n * 4, 0, sd), n, 4)
  list(X = X, Y = Y, z = z)
}

test_that("RV coefficient is high for shared, low for independent, with a valid p", {
  d <- shared_XY()
  rv <- rvCoefficient(d$X, d$Y, permutations = 500, seed = 1)
  expect_s3_class(rv, "rv_coefficient")
  expect_gt(rv$rv, 0.5); expect_true(rv$rv >= 0 && rv$rv <= 1)
  expect_lt(rv$p_value, 0.01)
  set.seed(2); Xi <- matrix(rnorm(50 * 5), 50, 5); Yi <- matrix(rnorm(50 * 4), 50, 4)
  rv0 <- rvCoefficient(Xi, Yi, permutations = 500, seed = 2)
  expect_lt(rv0$rv, rv$rv); expect_gt(rv0$p_value, 0.05)
})

test_that("distance correlation detects NON-linear dependence", {
  set.seed(3); x <- rnorm(80); y <- x^2 + rnorm(80, 0, 0.2)   # quadratic: Pearson ~ 0
  dc <- distanceCorrelation(x, y, permutations = 500, seed = 3)
  expect_s3_class(dc, "distance_correlation")
  expect_gt(dc$dcor, 0.4)                              # dCor sees the dependence
  expect_lt(dc$p_value, 0.01)
  expect_lt(abs(cor(x, y)), 0.3)                       # ...that Pearson misses
})

test_that("distance correlation is ~0 for independent variables", {
  set.seed(4); x <- rnorm(100); y <- rnorm(100)
  dc <- distanceCorrelation(x, y, permutations = 500, seed = 4)
  expect_lt(dc$dcor, 0.3); expect_gt(dc$p_value, 0.05)
})

test_that("RSA finds matching observation-geometry across modalities", {
  d <- shared_XY(seed = 5)
  rsa <- representationalSimilarity(list(kin = d$X, emg = d$Y),
                                    permutations = 500, seed = 5)
  expect_s3_class(rsa, "rsa")
  expect_equal(dim(rsa$similarity), c(2, 2))
  expect_gt(rsa$similarity[1, 2], 0.4)                # shared geometry
  expect_lt(rsa$p_value, 0.01)
  # scrambling one block's rows destroys the representational similarity
  set.seed(6); rsa0 <- representationalSimilarity(list(a = d$X, b = d$Y[sample(nrow(d$Y)), ]))
  expect_lt(rsa0$similarity[1, 2], rsa$similarity[1, 2])
})

test_that("RSA handles > 2 blocks and the correlation RDM", {
  d <- shared_XY(seed = 7)
  W <- d$z %*% matrix(rnorm(2 * 3), 2, 3) + matrix(rnorm(nrow(d$z) * 3, 0, .3), nrow(d$z), 3)
  rsa <- representationalSimilarity(list(a = d$X, b = d$Y, c = W), method = "correlation")
  expect_equal(dim(rsa$similarity), c(3, 3))
  expect_true(all(diag(rsa$similarity) == 1))
  expect_output(print(rsa), "Representational similarity")
})

test_that("mismatched rows error", {
  expect_error(rvCoefficient(matrix(0, 10, 3), matrix(0, 8, 3)), "same number of rows")
  expect_error(distanceCorrelation(matrix(0, 10, 3), matrix(0, 8, 3)), "same number of rows")
})
