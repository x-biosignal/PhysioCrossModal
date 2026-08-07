test_that("Wilson plus operator keeps causal lags and triangular zero lag", {
  set.seed(811)
  n_fft <- 32L
  lag <- array(
    complex(real = rnorm(2 * 2 * n_fft), imaginary = rnorm(2 * 2 * n_fft)),
    dim = c(2L, 2L, n_fft)
  )
  frequency <- PhysioCrossModal:::.fft_matrix_array(lag)
  projected <- PhysioCrossModal:::.wilson_plus(frequency)
  got <- PhysioCrossModal:::.fft_matrix_array(projected, inverse = TRUE)

  expected_zero <- lag[, , 1L]
  expected_zero[lower.tri(expected_zero)] <- 0
  diag(expected_zero) <- diag(expected_zero) / 2
  expect_equal(got[, , 1L], expected_zero, tolerance = 1e-12)
  expect_equal(
    got[, , 2:(n_fft %/% 2L)],
    lag[, , 2:(n_fft %/% 2L)],
    tolerance = 1e-12
  )
  expect_lt(max(Mod(got[, , (n_fft %/% 2L + 1L):n_fft])), 1e-12)
})

test_that("Wilson factorization recovers an analytic minimum-phase factor", {
  n_fft <- 128L
  a0 <- matrix(c(1, 0, 0.2, 0.8), 2L, 2L)
  a1 <- matrix(c(0.2, 0.1, 0, 0.1), 2L, 2L)
  psi <- array(0 + 0i, dim = c(2L, 2L, n_fft))
  for (frequency in 0:(n_fft - 1L)) {
    psi[, , frequency + 1L] <-
      a0 + a1 * exp(-2i * pi * frequency / n_fft)
  }
  spectrum <- array(0 + 0i, dim = dim(psi))
  expected_h <- array(0 + 0i, dim = dim(psi))
  for (frequency in seq_len(n_fft)) {
    spectrum[, , frequency] <- psi[, , frequency] %*%
      Conj(t(psi[, , frequency]))
    expected_h[, , frequency] <- psi[, , frequency] %*% solve(a0)
  }

  fit <- PhysioCrossModal:::.wilson_sfactor(
    spectrum, tolerance = 1e-10, max_iterations = 100L
  )
  expect_true(fit$converged)
  expect_lt(fit$reconstruction_error, 1e-12)
  expect_equal(fit$Sigma, a0 %*% t(a0), tolerance = 1e-12)
  expect_equal(fit$H, expected_h, tolerance = 1e-11)

  non_hermitian <- spectrum
  non_hermitian[1L, 2L, 1L] <- non_hermitian[1L, 2L, 1L] + 1i
  expect_error(
    PhysioCrossModal:::.wilson_sfactor(
      non_hermitian, tolerance = 0.9, max_iterations = 10L
    ),
    "not Hermitian"
  )
})

test_that("nonparametric Granger detects direction and matches magnitude", {
  set.seed(42)
  signals <- make_directed_signals(
    n = 5000, sr = 500, lag_samples = 5, coupling = 0.6
  )
  nonparametric <- grangerCausality(
    signals$x, signals$y, sr = signals$sr, method = "nonparametric"
  )
  parametric <- grangerCausality(
    signals$x, signals$y, sr = signals$sr, order = 10L
  )

  expect_gt(nonparametric$gc_xy, nonparametric$gc_yx)
  expect_gt(nonparametric$net_gc, 0)
  expect_lt(
    abs(nonparametric$gc_xy - parametric$gc_xy) / parametric$gc_xy,
    0.2
  )
  expect_identical(nonparametric$order, NA_integer_)
  expect_identical(nonparametric$method, "nonparametric")
  expect_true(nonparametric$converged)
  expect_lt(nonparametric$factorization_error, 0.02)
  expect_lte(
    nonparametric$factorization_error,
    nonparametric$factorization_limit
  )
  expect_lt(nonparametric$regularization_adjustment, 1e-8)
  expect_equal(
    length(nonparametric$frequencies),
    length(nonparametric$gc_xy_spectrum)
  )
  expect_true(all(nonparametric$gc_xy_spectrum >= 0))
  expect_true(all(nonparametric$gc_yx_spectrum >= 0))

  independent_xy <- PhysioCrossModal:::.gc_trapezoid_mean(
    nonparametric$gc_xy_spectrum,
    nonparametric$frequencies,
    nonparametric$freq_range
  )
  expect_equal(nonparametric$gc_xy, independent_xy, tolerance = 1e-12)

  band <- grangerCausality(
    signals$x, signals$y, sr = signals$sr, method = "nonparametric",
    freq_range = c(10, 40)
  )
  expect_equal(
    band$gc_xy,
    PhysioCrossModal:::.gc_trapezoid_mean(
      band$gc_xy_spectrum, band$frequencies, c(10, 40)
    ),
    tolerance = 1e-12
  )
})

test_that("directional spectra swap and survive scale or sign changes", {
  set.seed(410)
  signals <- make_directed_signals(
    n = 3000, sr = 500, lag_samples = 5, coupling = 0.6
  )
  base <- grangerCausality(
    signals$x, signals$y, sr = signals$sr, method = "nonparametric"
  )
  swapped <- grangerCausality(
    signals$y, signals$x, sr = signals$sr, method = "nonparametric"
  )
  transformed <- grangerCausality(
    7 * signals$x, -0.3 * signals$y,
    sr = signals$sr, method = "nonparametric"
  )

  expect_equal(base$gc_xy, swapped$gc_yx, tolerance = 1e-5)
  expect_equal(base$gc_yx, swapped$gc_xy, tolerance = 1e-5)
  expect_equal(base$gc_xy_spectrum, swapped$gc_yx_spectrum, tolerance = 2e-4)
  expect_equal(base$gc_yx_spectrum, swapped$gc_xy_spectrum, tolerance = 2e-4)
  expect_equal(base$gc_xy, transformed$gc_xy, tolerance = 1e-5)
  expect_equal(base$gc_yx, transformed$gc_yx, tolerance = 1e-5)
})

test_that("Welch and multitaper estimators preserve simulated direction", {
  set.seed(71)
  signals <- make_directed_signals(
    n = 5000, sr = 500, lag_samples = 5, coupling = 0.6
  )
  multitaper <- grangerCausality(
    signals$x, signals$y, sr = signals$sr, method = "nonparametric"
  )
  welch <- grangerCausality(
    signals$x, signals$y, sr = signals$sr, method = "nonparametric",
    spectral_estimator = "welch"
  )

  expect_gt(multitaper$gc_xy, multitaper$gc_yx)
  expect_gt(welch$gc_xy, welch$gc_yx)
  expect_gte(multitaper$n_realizations, 2L)
  expect_gte(welch$n_realizations, 2L)
  expect_identical(welch$spectral_estimator, "welch")
  expect_lt(welch$factorization_error, 0.02)
})

test_that("nonparametric Granger matches the pinned FieldTrip reference", {
  fixture <- readRDS(test_path(
    "fixtures", "nonparametric_granger_fieldtrip_reference.rds"
  ))
  estimated <- PhysioCrossModal:::.estimate_gc_csd(
    fixture$x - mean(fixture$x),
    fixture$y - mean(fixture$y),
    fixture$sr,
    "multitaper",
    fixture$n_fft,
    fixture$segment_length,
    fixture$overlap,
    fixture$time_bandwidth,
    fixture$n_tapers
  )
  reference_csd <- fixture$csd_real + 1i * fixture$csd_imag
  csd_nrmse <- sqrt(mean(Mod(estimated$spectrum - reference_csd)^2)) /
    sqrt(mean(Mod(reference_csd)^2))
  expect_lt(csd_nrmse, 1e-10)
  expect_identical(estimated$n_realizations, fixture$n_realizations)

  got <- grangerCausality(
    fixture$x,
    fixture$y,
    sr = fixture$sr,
    method = "nonparametric",
    n_fft = fixture$n_fft,
    segment_length = fixture$segment_length,
    overlap = fixture$overlap,
    time_bandwidth = fixture$time_bandwidth,
    n_tapers = fixture$n_tapers
  )
  scale_xy <- max(diff(range(fixture$gc_xy)), sqrt(.Machine$double.eps))
  scale_yx <- max(diff(range(fixture$gc_yx)), sqrt(.Machine$double.eps))
  nrmse_xy <- sqrt(mean((got$gc_xy_spectrum - fixture$gc_xy)^2)) / scale_xy
  nrmse_yx <- sqrt(mean((got$gc_yx_spectrum - fixture$gc_yx)^2)) / scale_yx
  expect_lt(nrmse_xy, 0.05)
  expect_lt(nrmse_yx, 0.05)
  expect_equal(
    got$gc_xy,
    PhysioCrossModal:::.gc_trapezoid_mean(
      fixture$gc_xy, fixture$frequencies, c(0, fixture$sr / 2)
    ),
    tolerance = 0.005
  )
  expect_equal(
    got$gc_yx,
    PhysioCrossModal:::.gc_trapezoid_mean(
      fixture$gc_yx, fixture$frequencies, c(0, fixture$sr / 2)
    ),
    tolerance = 0.005
  )
  expect_identical(
    fixture$provenance$fieldtrip_revision,
    "4c553bda8fb238b0afc549c90388d9c86124d126"
  )
  expect_identical(fixture$provenance$octave_version, "9.4.0")
})

test_that("independent signals have small nonparametric Granger values", {
  set.seed(993)
  result <- grangerCausality(
    rnorm(3000), rnorm(3000), sr = 500, method = "nonparametric"
  )
  expect_lt(result$gc_xy, 0.06)
  expect_lt(result$gc_yx, 0.06)
})

test_that("couplingAnalysis exposes a distinct Granger estimator selector", {
  set.seed(19)
  signals <- make_directed_signals(
    n = 3000, sr = 500, lag_samples = 5, coupling = 0.6
  )
  result <- couplingAnalysis(
    signals$x, signals$y, method = "granger", sr = signals$sr,
    granger_method = "nonparametric"
  )
  expect_identical(result$method, "nonparametric")
  expect_gt(result$gc_xy, result$gc_yx)

  parametric <- couplingAnalysis(
    signals$x, signals$y, method = "granger", sr = signals$sr
  )
  expect_identical(parametric$method, "parametric")
})

test_that("nonparametric Granger validates spectra and controls", {
  set.seed(210)
  x <- rnorm(1000)
  y <- rnorm(1000)
  call_gc <- function(...) {
    grangerCausality(x, y, sr = 200, method = "nonparametric", ...)
  }

  expect_error(call_gc(n_fft = 255L), "n_fft")
  expect_error(call_gc(segment_length = 2000L), "segment_length")
  expect_error(call_gc(overlap = 1), "overlap")
  expect_error(call_gc(time_bandwidth = 1), "time_bandwidth")
  expect_error(call_gc(n_tapers = 1L), "n_tapers")
  expect_error(call_gc(factor_tolerance = 0), "factor_tolerance")
  expect_error(call_gc(factor_max_iterations = 1L), "did not converge")
  expect_error(
    call_gc(factor_tolerance = 0.9),
    "reconstruction error"
  )
  expect_error(call_gc(eigen_floor = 0), "eigen_floor")
  expect_error(call_gc(freq_range = c(-1, 20)), "freq_range")
  expect_error(call_gc(freq_range = c(20, 101)), "freq_range")
  expect_error(call_gc(freq_range = c(1.1, 1.4)), "two spectral bins")
  expect_error(
    grangerCausality(x, y, sr = 0, method = "nonparametric"),
    "sampling rate|sr"
  )
  expect_error(
    grangerCausality(x[seq_len(31)], y[seq_len(31)], sr = 200,
                     method = "nonparametric"),
    "at least 32"
  )
  expect_error(
    grangerCausality(x[seq_len(100)], y[seq_len(100)], sr = 200,
                     method = "nonparametric", spectral_estimator = "welch",
                     segment_length = 80L, overlap = 0),
    "two complete segments"
  )
  expect_error(
    grangerCausality(x, y[-1], sr = 200, method = "nonparametric"),
    "equal-length"
  )
  x_bad <- x
  x_bad[3] <- NA_real_
  expect_error(
    grangerCausality(x_bad, y, sr = 200, method = "nonparametric"),
    "finite"
  )
  expect_error(
    grangerCausality(x, x, sr = 200, method = "nonparametric"),
    "rank-deficient"
  )
  expect_error(
    grangerCausality(rep(1, 1000), y, sr = 200, method = "nonparametric"),
    "non-constant"
  )
})

test_that("spectral regularization is bounded and rejects indefiniteness", {
  spectrum <- array(0 + 0i, dim = c(2L, 2L, 3L))
  for (frequency in seq_len(3L)) {
    spectrum[, , frequency] <- diag(c(1, 1e-14))
  }
  adjusted <- PhysioCrossModal:::.regularize_gc_spectrum(spectrum, 1e-10)
  expect_gt(adjusted$maximum_adjustment, 0)
  expect_lt(adjusted$maximum_adjustment, 1e-8)

  spectrum[, , 2L] <- diag(c(1, -1e-4))
  expect_error(
    PhysioCrossModal:::.regularize_gc_spectrum(spectrum, 1e-10),
    "materially indefinite"
  )
})

test_that("trapezoidal aggregation uses endpoint half weights", {
  frequencies <- 0:4
  values <- frequencies^2
  got <- PhysioCrossModal:::.gc_trapezoid_mean(
    values, frequencies, c(0, 4)
  )
  expected <- sum(diff(frequencies) *
                    (head(values, -1L) + tail(values, -1L)) / 2) / 4
  expect_equal(got, expected, tolerance = 1e-15)
  expect_error(
    PhysioCrossModal:::.gc_trapezoid_mean(values, frequencies, c(1.1, 1.9)),
    "two spectral bins"
  )
})
