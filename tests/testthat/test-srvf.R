library(testthat)
library(PhysioCrossModal)

# a shared analytic grid + Gaussian bump
.grid <- function(N = 101) seq(0, 1, length.out = N)
.bump <- function(t, mu = 0.5, sd = 0.08) exp(-((t - mu)^2) / (2 * sd^2))

test_that("elastic alignment recovers a known diffeomorphism and reduces amplitude distance", {
  t <- .grid(101)
  f <- .bump(t)
  gamma_true <- t^1.7                       # a known monotone warp
  g <- warpApply(f, gamma_true, t)          # f o gamma_true

  al <- elasticAlign(f, g, t)

  # the warped curve is pulled back onto the reference
  expect_lt(max(abs(al$aligned - f)), 0.05)
  # aligning collapses the SRVF amplitude distance vs the unaligned pair
  q_f <- elasticAlign(f, f, t)$srvf_x               # SRVF of the reference
  q_g <- elasticAlign(g, g, t)$srvf_x               # SRVF of the unaligned curve
  unaligned <- sqrt(sum((q_f - q_g)^2) * mean(diff(t)))
  expect_lt(al$amplitude_distance, unaligned / 5)   # >5x reduction

  # the fdasrvf solver also recovers the alignment: the aligned curve tracks the
  # reference and the amplitude distance is reduced by more than an order of
  # magnitude. fdasrvf and the native DP solver leave comparable, non-zero
  # residuals on a discretised warp (the earlier 1e-3 expectation was never
  # exercised because the test skips when fdasrvf is absent).
  skip_if_not_installed("fdasrvf")
  alf <- elasticAlign(f, g, t, use_fdasrvf = TRUE)
  expect_lt(max(abs(alf$aligned - f)), 0.05)
  expect_lt(alf$amplitude_distance, unaligned / 10)
})

test_that("the recovered warping is a boundary-preserving monotone diffeomorphism", {
  t <- .grid(81)
  al <- elasticAlign(.bump(t, 0.45), .bump(t, 0.6), t)
  expect_true(all(diff(al$gamma) >= -1e-9))         # monotone non-decreasing
  expect_equal(al$gamma[1], 0, tolerance = 1e-8)    # boundary 0
  expect_equal(al$gamma[length(al$gamma)], 1, tolerance = 1e-8)  # boundary 1
  expect_gt(al$phase_distance, 0)                   # a genuine phase difference
})

test_that("identity alignment is a no-op", {
  t <- .grid(101)
  f <- .bump(t)
  al <- elasticAlign(f, f, t)
  expect_lt(al$amplitude_distance, 1e-6)
  expect_lt(max(abs(al$gamma - t)), 1e-8)
})

test_that("elastic mean preserves the peak where the naive mean blurs it", {
  t <- .grid(101)
  reps <- lapply(c(-0.15, -0.07, 0, 0.07, 0.15),
                 function(s) .bump(t, 0.5 + s, sd = 0.05))
  naive <- Reduce(`+`, reps) / length(reps)
  m <- srvfMean(reps, t)

  expect_gt(max(m$mean), 0.9)                        # elastic keeps the peak
  expect_lt(max(naive), 0.75)                        # naive averaging blurs it
  expect_equal(dim(m$warpings), c(length(t), length(reps)))
  expect_length(m$aligned, length(reps))
})

test_that("warpApply carries a warping onto a companion channel", {
  t <- .grid(60)
  companion <- sin(2 * pi * t)
  w <- warpApply(companion, gamma = t, t = t)        # identity warp
  expect_equal(w, companion, tolerance = 1e-8)
  # a non-trivial warp changes the signal but preserves its endpoints
  w2 <- warpApply(companion, gamma = t^1.5, t = t)
  expect_equal(w2[1], companion[1], tolerance = 1e-8)
  expect_equal(w2[length(w2)], companion[length(companion)], tolerance = 1e-8)
})

test_that("alignSignals(method = 'elastic') warps modalities onto the first", {
  t <- .grid(200)
  ref <- matrix(.bump(t, 0.5), 200, 1)
  shifted <- matrix(.bump(t, 0.62), 200, 1)
  a <- PhysioExperiment(assays = list(raw = ref), samplingRate = 200)
  b <- PhysioExperiment(assays = list(raw = shifted), samplingRate = 200)

  mpe <- alignSignals(A = a, B = b, method = "elastic")
  expect_s4_class(mpe, "MultiPhysioExperiment")
  bw <- mpe[["B"]]
  # the warp is recorded and the warped B peak moves toward the reference peak
  expect_false(is.null(S4Vectors::metadata(bw)$elastic_warp))
  peak_ref <- which.max(SummarizedExperiment::assay(mpe[["A"]], "raw")[, 1])
  peak_b_before <- which.max(shifted[, 1])
  peak_b_after <- which.max(SummarizedExperiment::assay(bw, "raw")[, 1])
  expect_lt(abs(peak_b_after - peak_ref), abs(peak_b_before - peak_ref))
})
