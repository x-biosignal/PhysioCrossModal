library(testthat)
library(PhysioCrossModal)

# Gamma (fa) amplitude coupled to the phase of a narrowband (jittery) theta
# oscillation at fp; coupling = 0 is uncoupled. A jittery phase (rather than a
# pure sine) is the realistic case where timeshift surrogates break coupling.
make_pac_signal <- function(n = 5000, sr = 500, fp = 6, fa = 50, coupling = 1,
                            seed = 1) {
  set.seed(seed)
  t <- (0:(n - 1)) / sr
  theta <- suppressWarnings(
    PhysioCrossModal:::.bandpass_filter(stats::rnorm(n), fp - 1, fp + 1, sr))
  theta <- theta / stats::sd(theta)
  ph <- PhysioCrossModal:::.hilbert_phase(theta)
  modu <- (1 + cos(ph - pi)) / 2                 # in 0..1
  theta + (0.1 + coupling * modu) * sin(2 * pi * fa * t) + stats::rnorm(n, sd = 0.3)
}

test_that("Tort MI is high for strong coupling and near 0 when uncoupled", {
  sc <- make_pac_signal(coupling = 1)
  su <- make_pac_signal(coupling = 0)
  mi_c <- phaseAmplitudeCoupling(sc, sr = 500, phase_band = c(4, 8),
                                 amp_band = c(35, 65), method = "tort")$pac
  mi_u <- phaseAmplitudeCoupling(su, sr = 500, phase_band = c(4, 8),
                                 amp_band = c(35, 65), method = "tort")$pac
  expect_gt(mi_c, 0.01)
  expect_lt(mi_u, 0.005)
  expect_gt(mi_c, 10 * mi_u)
})

test_that("Tort MI is in 0..1 and ~0 for a uniform phase-amplitude distribution", {
  set.seed(2)
  phase <- stats::runif(5000, -pi, pi)
  amp <- abs(stats::rnorm(5000)) + 1        # amplitude independent of phase
  mi <- PhysioCrossModal:::.pac_tort(phase, amp)$mi
  expect_gte(mi, 0)
  expect_lte(mi, 1)
  expect_lt(mi, 0.005)
})

test_that("Canolty and Ozkurt measures track the coupling; Ozkurt is bounded", {
  sc <- make_pac_signal(coupling = 1); su <- make_pac_signal(coupling = 0)
  oz_c <- phaseAmplitudeCoupling(sc, sr = 500, phase_band = c(4, 8),
                                 amp_band = c(35, 65), method = "ozkurt")$pac
  oz_u <- phaseAmplitudeCoupling(su, sr = 500, phase_band = c(4, 8),
                                 amp_band = c(35, 65), method = "ozkurt")$pac
  expect_gt(oz_c, oz_u); expect_gte(oz_u, 0); expect_lte(oz_c, 1)
  ca_c <- phaseAmplitudeCoupling(sc, sr = 500, phase_band = c(4, 8),
                                 amp_band = c(35, 65), method = "canolty")$pac
  expect_gt(ca_c, 0)
})

test_that("comodulogram peaks at the true (phase, amplitude) frequencies", {
  sig <- make_pac_signal(fp = 6, fa = 50, coupling = 1)
  cm <- comodulogram(sig, sr = 500, phase_freqs = c(2, 4, 6, 8, 10),
                     amp_freqs = seq(20, 80, by = 10), method = "tort")
  expect_equal(dim(cm$matrix), c(5L, 7L))
  expect_equal(cm$peak$phase_freq, 6)
  expect_equal(cm$peak$amp_freq, 50)
})

test_that("plotComodulogram returns a ggplot", {
  skip_if_not_installed("ggplot2")
  sig <- make_pac_signal()
  cm <- comodulogram(sig, sr = 500, phase_freqs = c(4, 6, 8),
                     amp_freqs = c(30, 50, 70))
  expect_s3_class(plotComodulogram(cm), "ggplot")
})

test_that("surrogate test: significant for coupled, non-significant for uncoupled", {
  sc <- make_pac_signal(coupling = 1)
  su <- make_pac_signal(coupling = 0)
  p_c <- surrogateTest(sc, sc, sr = 500, method = "pac", n_surrogates = 99,
                       surrogate_type = "timeshift", phase_band = c(4, 8),
                       amp_band = c(35, 65))$p_value
  p_u <- surrogateTest(su, su, sr = 500, method = "pac", n_surrogates = 99,
                       surrogate_type = "timeshift", phase_band = c(4, 8),
                       amp_band = c(35, 65))$p_value
  expect_lt(p_c, 0.05)
  expect_gt(p_u, 0.05)
})

test_that("phaseAmplitudeCoupling validates bands and supports two-signal input", {
  sig <- make_pac_signal()
  expect_error(phaseAmplitudeCoupling(sig, sr = 500, phase_band = c(8, 4)),
               "phase_band")
  # cross-signal: phase from x, amplitude from y
  r <- phaseAmplitudeCoupling(sig, sig, sr = 500, phase_band = c(4, 8),
                              amp_band = c(35, 65))
  expect_true(is.finite(r$pac))
})
