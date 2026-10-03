# Composite multi-modal gait index, feature fusion, cross-block factor model.

test_that("multimodal gait index scores controls ~100 and deviants lower", {
  set.seed(1); p <- 8
  ctrl <- matrix(rnorm(80 * p), 80, p)
  typical <- matrix(rnorm(6 * p), 6, p)
  deviant <- matrix(rnorm(6 * p, 4), 6, p)             # far from the control cloud
  idx <- multimodalGaitIndex(rbind(typical, deviant), ctrl)
  expect_s3_class(idx, "multimodal_index")
  expect_equal(idx$reference_mean_index, 100, tolerance = 2)
  expect_gt(mean(idx$index[1:6]), mean(idx$index[7:12]))    # typical > deviant
  expect_lt(mean(idx$index[7:12]), 90)                      # deviants clearly below 100
})

test_that("multimodal gait index breaks the deviation down by modality", {
  set.seed(2); n_ctrl <- 100
  ctrl <- matrix(rnorm(n_ctrl * 9), n_ctrl, 9)
  blocks <- list(kinematics = 1:3, kinetics = 4:6, emg = 7:9)
  # a subject deviant ONLY in the EMG block
  subj <- matrix(rnorm(9), 1, 9); subj[7:9] <- subj[7:9] + 5
  idx <- multimodalGaitIndex(subj, ctrl, blocks = blocks)
  expect_true(!is.null(idx$by_modality))
  expect_equal(colnames(idx$by_modality), c("kinematics", "kinetics", "emg"))
  # EMG modality index is the lowest (most deviant)
  expect_equal(which.min(idx$by_modality[1, ]), c(emg = 3L), ignore_attr = TRUE)
})

test_that("fuseBlocks balances blocks of very different scale", {
  set.seed(3)
  small <- matrix(rnorm(20 * 3), 20, 3)
  big   <- matrix(rnorm(20 * 2) * 1000, 20, 2)
  f <- fuseBlocks(list(small = small, big = big), method = "zscore")
  expect_s3_class(f, "fused_blocks")
  expect_equal(dim(f$matrix), c(20, 5))
  # after z-scoring, every fused column has unit-ish variance (no block dominates)
  expect_true(all(abs(apply(f$matrix, 2, sd) - 1) < 0.01))
  expect_equal(nrow(f$block_index), 2)
  expect_output(print(f), "Fused blocks")
})

test_that("fuseBlocks mfa/frobenius scalings run and preserve dimensions", {
  b <- list(a = matrix(rnorm(30), 10, 3), b = matrix(rnorm(20), 10, 2))
  for (m in c("mfa", "frobenius"))
    expect_equal(dim(fuseBlocks(b, method = m)$matrix), c(10, 5))
})

test_that("cross-block factor model flags a factor shared across modalities", {
  set.seed(4); n <- 120
  f <- matrix(rnorm(n), n, 1)                          # one common latent
  a <- f %*% matrix(rnorm(4), 1, 4) + matrix(rnorm(n * 4, 0, 0.4), n, 4)
  b <- f %*% matrix(rnorm(3), 1, 3) + matrix(rnorm(n * 3, 0, 0.4), n, 3)
  cf <- crossBlockFactor(list(kin = a, emg = b), nfactors = 1)
  expect_s3_class(cf, "crossblock_factor")
  expect_true(cf$shared[1])                            # factor 1 loads on both blocks
  expect_equal(cf$n_shared, 1)
  expect_output(print(cf), "Cross-block factor")
})

test_that("cross-block factor model marks a modality-specific factor as not shared", {
  set.seed(5); n <- 150
  fc <- matrix(rnorm(n), n, 1)                         # common factor
  fa <- matrix(rnorm(n), n, 1)                         # block-a-only factor
  a <- fc %*% matrix(rnorm(4), 1, 4) + fa %*% matrix(rnorm(4), 1, 4) +
       matrix(rnorm(n * 4, 0, 0.3), n, 4)
  b <- fc %*% matrix(rnorm(4), 1, 4) + matrix(rnorm(n * 4, 0, 0.3), n, 4)
  cf <- crossBlockFactor(list(a = a, b = b), nfactors = 2)
  expect_true(any(cf$shared))                          # the common factor is shared
  expect_true(any(!cf$shared))                         # the block-a factor is not
})

test_that("index and fusion input checks error", {
  expect_error(multimodalGaitIndex(matrix(0, 5, 3), matrix(0, 5, 4)), "same columns")
  expect_error(fuseBlocks(list(matrix(0, 5, 2))), ">= 2")
})
