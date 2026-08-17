library(testthat)
library(PhysioCrossModal)

# 2-channel PE with a unidirectional x -> y coupling (1-sample lag).
make_chain2 <- function(n = 2000, sr = 100, seed = 1) {
  set.seed(seed)
  x <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
  y <- numeric(n); ny <- stats::rnorm(n, sd = 0.4)
  for (t in 2:n) y[t] <- 0.6 * x[t - 1] + 0.2 * y[t - 1] + ny[t]
  PhysioExperiment(
    assays = list(raw = cbind(x, y)),
    colData = S4Vectors::DataFrame(label = c("x", "y"), type = rep("sig", 2)),
    samplingRate = sr)
}

# Unidirectional source -> target with a 5-sample lag.
lagged_pair <- function(n = 1500, seed = 1) {
  set.seed(seed)
  x <- stats::rnorm(n)
  y <- c(rep(0, 5), 0.7 * x[1:(n - 5)]) + stats::rnorm(n, sd = 0.5)
  list(x = x, y = y)
}

test_that("transfer entropy is directional for unidirectional x -> y (5-sample lag)", {
  nets <- vapply(1:4, function(s) {
    d <- lagged_pair(seed = s)
    transferEntropy(d$x, d$y, sr = 100, delay = 5, estimator = "histogram")$net
  }, numeric(1))
  expect_true(all(nets > 0))                        # x->y dominates every seed
  d <- lagged_pair(seed = 1)
  r <- transferEntropy(d$x, d$y, sr = 100, delay = 5, estimator = "histogram")
  expect_gt(r$te_xy, 3 * abs(r$te_yx))
})

test_that("KSG transfer entropy also recovers the direction", {
  d <- lagged_pair(n = 700, seed = 1)
  r <- transferEntropy(d$x, d$y, sr = 100, delay = 5, estimator = "ksg", knn = 4)
  expect_gt(r$te_xy, r$te_yx)
  expect_gt(r$net, 0)
})

test_that("effective transfer entropy of independent signals is within the null band", {
  set.seed(3); a <- stats::rnorm(1500); b <- stats::rnorm(1500)
  r <- transferEntropy(a, b, sr = 100, delay = 1, estimator = "histogram",
                       effective = TRUE, n_surrogate = 15, seed = 7)
  expect_lt(abs(r$eff_xy), 3 * r$null_sd_xy)
  expect_lt(abs(r$eff_yx), 3 * r$null_sd_yx)
})

test_that("DTF / PDC match PhysioEEG eegDTF / eegPDC (parity < 1e-6)", {
  skip_if_not_installed("PhysioEEG")
  if (!exists("eegDTF", where = asNamespace("PhysioEEG"))) skip("eegDTF unavailable")
  pe <- make_chain2()
  arr_cm <- directedTransferFunction(pe, order = 5)$dtf
  arr_eeg <- S4Vectors::metadata(PhysioEEG::eegDTF(pe, order = 5))$connectivity$array
  expect_lt(max(abs(arr_cm - arr_eeg)), 1e-6)

  parr_cm <- partialDirectedCoherence(pe, order = 5)$pdc
  parr_eeg <- S4Vectors::metadata(PhysioEEG::eegPDC(pe, order = 5))$connectivity$array
  expect_lt(max(abs(parr_cm - parr_eeg)), 1e-6)
})

test_that("DTF captures the indirect path where PDC does not (3-node chain)", {
  set.seed(1); n <- 3000; e <- matrix(stats::rnorm(n * 3), n, 3); X <- matrix(0, n, 3)
  for (t in 2:n) {
    X[t, 1] <- 0.4 * X[t - 1, 1] + e[t, 1]
    X[t, 2] <- 0.5 * X[t - 1, 1] + e[t, 2]
    X[t, 3] <- 0.5 * X[t - 1, 2] + e[t, 3]
  }
  d <- directedTransferFunction(X, sr = 100, order = 4)$dtf
  p <- partialDirectedCoherence(X, sr = 100, order = 4)$pdc
  expect_gt(mean(d[3, 1, ]), 0.1)                   # DTF 1->3 (indirect)
  expect_lt(mean(p[3, 1, ]), 0.1)                   # PDC 1->3 (~0)
})

test_that("couplingMatrix dispatches dtf, pdc, and transferentropy", {
  pe <- make_chain2()
  m_dtf <- couplingMatrix(pe, pe, method = "dtf", order = 5)$matrix
  expect_true(is.matrix(m_dtf))
  expect_true(is.finite(m_dtf["x", "y"]))        # off-diagonal x->y computed
  expect_true(is.matrix(couplingMatrix(pe, pe, method = "pdc", order = 5)$matrix))
  m_te <- couplingMatrix(pe, pe, method = "transferentropy",
                         estimator = "histogram", delay = 1)$matrix
  expect_true(is.matrix(m_te))
  expect_true(is.finite(m_te["x", "y"]))
})

test_that("directedTransferFunction validates input", {
  expect_error(directedTransferFunction(matrix(stats::rnorm(300), ncol = 3)), "sr")
  expect_error(directedTransferFunction(stats::rnorm(100)), "Provide")
})
