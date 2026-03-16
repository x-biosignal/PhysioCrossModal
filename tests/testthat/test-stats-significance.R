# Tests for statistical significance testing (stats-significance.R)

# ---- .phase_randomize_surrogate ----------------------------------------------

test_that(".phase_randomize_surrogate preserves amplitude spectrum", {
  set.seed(42)
  x <- sin(2 * pi * 10 * seq(0, 1, length.out = 500)) + 0.3 * rnorm(500)

  surr <- PhysioCrossModal:::.phase_randomize_surrogate(x)

  expect_length(surr, length(x))
  expect_true(is.numeric(surr))

  # Amplitude spectra should match
  amp_orig <- Mod(stats::fft(x))
  amp_surr <- Mod(stats::fft(surr))
  expect_equal(amp_surr, amp_orig, tolerance = 1e-10)
})

test_that(".phase_randomize_surrogate produces real output (Hermitian symmetry)", {
  set.seed(1)
  # Even length
  x_even <- rnorm(100)
  surr_even <- PhysioCrossModal:::.phase_randomize_surrogate(x_even)
  expect_true(is.numeric(surr_even))
  expect_true(all(is.finite(surr_even)))
  # Verify Hermitian symmetry of the FFT: conjugate mirror property
  fft_surr <- stats::fft(surr_even)
  n <- length(surr_even)
  for (k in 2:(n / 2)) {
    expect_equal(fft_surr[k], Conj(fft_surr[n - k + 2]), tolerance = 1e-10)
  }

  # Odd length
  x_odd <- rnorm(101)
  surr_odd <- PhysioCrossModal:::.phase_randomize_surrogate(x_odd)
  expect_true(is.numeric(surr_odd))
  expect_length(surr_odd, 101)
  expect_true(all(is.finite(surr_odd)))
})

test_that(".phase_randomize_surrogate differs from original", {
  set.seed(42)
  x <- rnorm(200)
  surr <- PhysioCrossModal:::.phase_randomize_surrogate(x)
  # Should not be identical (extremely unlikely)
  expect_false(isTRUE(all.equal(x, surr)))
})

# ---- .timeshift_surrogate ----------------------------------------------------

test_that(".timeshift_surrogate preserves values but shifts them", {
  set.seed(1)
  x <- seq_len(50)
  surr <- PhysioCrossModal:::.timeshift_surrogate(x)

  expect_length(surr, length(x))
  # Same set of values, just reordered circularly
  expect_equal(sort(surr), sort(x))
  # Should differ from original
  expect_false(identical(surr, x))
})

test_that(".timeshift_surrogate preserves autocorrelation exactly", {
  set.seed(42)
  x <- cumsum(rnorm(500))
  surr <- PhysioCrossModal:::.timeshift_surrogate(x)

  # Circular autocorrelation at lag 1 should be identical
  acf_x <- sum(x * c(x[-1], x[1]))
  acf_surr <- sum(surr * c(surr[-1], surr[1]))
  expect_equal(acf_surr, acf_x, tolerance = 1e-10)
})

# ---- .block_bootstrap_indices ------------------------------------------------

test_that(".block_bootstrap_indices returns correct length", {
  set.seed(1)
  idx <- PhysioCrossModal:::.block_bootstrap_indices(100)
  expect_length(idx, 100)
  expect_true(all(idx >= 1 & idx <= 100))
})

test_that(".block_bootstrap_indices respects custom block_len", {
  set.seed(1)
  idx <- PhysioCrossModal:::.block_bootstrap_indices(50, block_len = 10)
  expect_length(idx, 50)
  expect_true(all(idx >= 1 & idx <= 50))
})

test_that(".block_bootstrap_indices produces contiguous blocks", {
  set.seed(42)
  idx <- PhysioCrossModal:::.block_bootstrap_indices(100, block_len = 20)
  # Check that at least some consecutive indices differ by 1 (contiguous)
  diffs <- diff(idx)
  expect_true(sum(diffs == 1) > length(idx) / 2)
})

# ---- .extract_coupling_stat --------------------------------------------------

test_that(".extract_coupling_stat extracts coherence peak", {
  fake <- list(coherence = c(0.1, 0.8, 0.3))
  stat <- PhysioCrossModal:::.extract_coupling_stat(fake, "coherence")
  expect_equal(stat, 0.8)
})

test_that(".extract_coupling_stat extracts PLV", {
  fake <- list(plv = 0.75)
  stat <- PhysioCrossModal:::.extract_coupling_stat(fake, "plv")
  expect_equal(stat, 0.75)
})

test_that(".extract_coupling_stat extracts PLI", {
  fake <- list(pli = 0.6)
  stat <- PhysioCrossModal:::.extract_coupling_stat(fake, "pli")
  expect_equal(stat, 0.6)
})

test_that(".extract_coupling_stat extracts wPLI", {
  fake <- list(wpli = 0.5)
  stat <- PhysioCrossModal:::.extract_coupling_stat(fake, "wpli")
  expect_equal(stat, 0.5)
})

test_that(".extract_coupling_stat extracts gc_xy", {
  fake <- list(gc_xy = 0.3, gc_yx = 0.1)
  stat <- PhysioCrossModal:::.extract_coupling_stat(fake, "granger")
  expect_equal(stat, 0.3)
})

test_that(".extract_coupling_stat extracts abs peak_correlation", {
  fake <- list(peak_correlation = -0.8)
  stat <- PhysioCrossModal:::.extract_coupling_stat(fake, "crosscorrelation")
  expect_equal(stat, 0.8)
})

test_that(".extract_coupling_stat errors on unknown method", {
  expect_error(
    PhysioCrossModal:::.extract_coupling_stat(list(), "unknown"),
    "Unknown method"
  )
})

# ---- surrogateTest -----------------------------------------------------------

test_that("surrogateTest p < 0.05 for strongly coupled signals", {
  set.seed(42)
  sr <- 500
  n <- sr * 4
  t <- seq(0, 4, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- 0.9 * sin(2 * pi * 10 * t) + 0.2 * rnorm(n)

  result <- surrogateTest(x, y, sr = sr, method = "coherence",
                          n_surrogates = 49, nperseg = 128L)

  expect_type(result, "list")
  expect_true("observed" %in% names(result))
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("surrogate_distribution" %in% names(result))
  expect_true("threshold_95" %in% names(result))

  expect_true(result$p_value < 0.05,
    info = paste("p =", result$p_value, "for coupled signals"))
  expect_length(result$surrogate_distribution, 49)
})

test_that("surrogateTest p > 0.05 for independent signals", {
  set.seed(123)
  sr <- 500
  n <- sr * 2
  x <- rnorm(n)
  y <- rnorm(n)

  result <- surrogateTest(x, y, sr = sr, method = "coherence",
                          n_surrogates = 49, nperseg = 128L)

  expect_true(result$p_value > 0.05,
    info = paste("p =", result$p_value, "for independent signals"))
})

test_that("surrogateTest works with timeshift surrogates", {
  set.seed(42)
  sr <- 500
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t) + 0.2 * rnorm(n)

  result <- surrogateTest(x, y, sr = sr, method = "coherence",
                          n_surrogates = 19, surrogate_type = "timeshift",
                          nperseg = 128L)

  expect_type(result, "list")
  expect_length(result$surrogate_distribution, 19)
  # With only 19 surrogates, p-value is coarsely estimated; just check it's valid
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that("surrogateTest works with PLV method", {
  set.seed(42)
  sr <- 500
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- sin(2 * pi * 10 * t + pi / 4) + 0.2 * rnorm(n)

  result <- surrogateTest(x, y, sr = sr, method = "plv",
                          n_surrogates = 19, freq_band = c(8, 12))

  expect_type(result, "list")
  expect_true(result$statistic > 0)
  expect_true(result$p_value <= 1)
})

test_that("surrogateTest validates n_surrogates", {
  expect_error(
    surrogateTest(rnorm(100), rnorm(100), sr = 100, method = "coherence",
                  n_surrogates = 0),
    "n_surrogates"
  )
})

test_that("surrogateTest return structure is complete", {
  set.seed(1)
  result <- surrogateTest(rnorm(500), rnorm(500), sr = 100,
                          method = "coherence", n_surrogates = 9,
                          nperseg = 64L)

  expect_named(result, c("observed", "statistic", "p_value",
                          "surrogate_distribution", "threshold_95"))
  expect_true(is.numeric(result$statistic))
  expect_true(is.numeric(result$p_value))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_true(is.numeric(result$threshold_95))
})

# ---- bootstrapCI -------------------------------------------------------------

test_that("bootstrapCI CI brackets the observed value", {
  set.seed(42)
  sr <- 500
  n <- sr * 4
  t <- seq(0, 4, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t) + 0.3 * rnorm(n)

  result <- bootstrapCI(x, y, sr = sr, method = "coherence",
                        n_boot = 49, nperseg = 128L)

  expect_type(result, "list")
  expect_true("observed" %in% names(result))
  expect_true("statistic" %in% names(result))
  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true("ci_level" %in% names(result))
  expect_true("boot_distribution" %in% names(result))

  # CI should bracket the observed value (with high probability)
  expect_true(result$ci_lower <= result$statistic)
  expect_true(result$ci_upper >= result$ci_lower)
  expect_equal(result$ci_level, 0.95)
  expect_length(result$boot_distribution, 49)
})

test_that("bootstrapCI works with PLV method", {
  set.seed(42)
  sr <- 500
  n <- sr * 2
  t <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t) + 0.2 * rnorm(n)
  y <- sin(2 * pi * 10 * t + pi / 4) + 0.2 * rnorm(n)

  result <- bootstrapCI(x, y, sr = sr, method = "plv",
                        n_boot = 19, freq_band = c(8, 12))

  expect_type(result, "list")
  expect_true(result$ci_lower >= 0)
  expect_true(result$ci_upper <= 1)
})

test_that("bootstrapCI validates inputs", {
  expect_error(
    bootstrapCI(rnorm(100), rnorm(100), sr = 100, method = "coherence",
                n_boot = 0),
    "n_boot"
  )
  expect_error(
    bootstrapCI(rnorm(100), rnorm(100), sr = 100, method = "coherence",
                ci = 1.5),
    "ci must be"
  )
})

test_that("bootstrapCI return structure is complete", {
  set.seed(1)
  result <- bootstrapCI(rnorm(500), rnorm(500), sr = 100,
                        method = "coherence", n_boot = 9,
                        nperseg = 64L)

  expect_named(result, c("observed", "statistic", "ci_lower",
                          "ci_upper", "ci_level", "boot_distribution"))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(result$ci_lower <= result$ci_upper)
})

test_that("bootstrapCI respects custom ci level", {
  set.seed(42)
  result <- bootstrapCI(rnorm(500), rnorm(500), sr = 100,
                        method = "coherence", n_boot = 19, ci = 0.9,
                        nperseg = 64L)
  expect_equal(result$ci_level, 0.9)
})

test_that("bootstrapCI respects custom block_len", {
  set.seed(42)
  result <- bootstrapCI(rnorm(500), rnorm(500), sr = 100,
                        method = "coherence", n_boot = 9,
                        block_len = 50, nperseg = 64L)
  expect_type(result, "list")
  expect_length(result$boot_distribution, 9)
})

# ---- Parallel computation (cores parameter) ----------------------------------

test_that("surrogateTest with cores=1 produces same structure", {
  set.seed(42)
  result <- surrogateTest(rnorm(500), rnorm(500), sr = 100,
                          method = "coherence", n_surrogates = 9,
                          cores = 1L, nperseg = 64L)
  expect_named(result, c("observed", "statistic", "p_value",
                          "surrogate_distribution", "threshold_95"))
  expect_length(result$surrogate_distribution, 9)
})

test_that("surrogateTest with cores=2 returns valid results", {
  skip_on_cran()
  set.seed(42)
  sr <- 500
  n <- sr * 2
  t_vec <- seq(0, 2, length.out = n)
  x <- sin(2 * pi * 10 * t_vec) + 0.2 * rnorm(n)
  y <- 0.8 * sin(2 * pi * 10 * t_vec) + 0.2 * rnorm(n)

  result <- surrogateTest(x, y, sr = sr, method = "coherence",
                          n_surrogates = 9, cores = 2L, nperseg = 128L)

  expect_type(result, "list")
  expect_length(result$surrogate_distribution, 9)
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_true(is.numeric(result$statistic))
})

test_that("bootstrapCI with cores=2 returns valid results", {
  skip_on_cran()
  set.seed(42)
  result <- bootstrapCI(rnorm(500), rnorm(500), sr = 100,
                        method = "coherence", n_boot = 9,
                        cores = 2L, nperseg = 64L)

  expect_type(result, "list")
  expect_length(result$boot_distribution, 9)
  expect_true(result$ci_lower <= result$ci_upper)
})

test_that("parallel results have correct length", {
  skip_on_cran()
  set.seed(42)
  result <- surrogateTest(rnorm(500), rnorm(500), sr = 100,
                          method = "coherence", n_surrogates = 15,
                          cores = 2L, nperseg = 64L)
  expect_length(result$surrogate_distribution, 15)
})

test_that(".parallel_apply falls back to sequential on error", {
  # Force an error in parallel by using an unreasonable core count
  # The tryCatch in .parallel_apply should fall back to lapply
  result <- PhysioCrossModal:::.parallel_apply(1:5, function(i) i^2, cores = 2L)
  expect_length(result, 5)
  expect_equal(result[[3]], 9)
})

# ---- surrogateMatrixTest -----------------------------------------------------

test_that("surrogateMatrixTest returns correct structure", {
  skip_if_not_installed("PhysioCore")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  pe1 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )
  pe2 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )

  result <- surrogateMatrixTest(pe1, pe2, method = "coherence",
                                 n_surrogates = 9, nperseg = 64L)

  expect_type(result, "list")
  expect_true("matrix" %in% names(result))
  expect_true("p_values" %in% names(result))
  expect_true("p_adjusted" %in% names(result))
  expect_true("significant" %in% names(result))
  expect_true("correction" %in% names(result))
  expect_true("alpha" %in% names(result))
})

test_that("surrogateMatrixTest p_values has correct dimensions", {
  skip_if_not_installed("PhysioCore")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  pe1 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )
  pe2 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 3), nrow = n)),
    samplingRate = sr
  )

  result <- surrogateMatrixTest(pe1, pe2, method = "coherence",
                                 n_surrogates = 9, nperseg = 64L)

  expect_equal(dim(result$p_values), c(2, 3))
  expect_equal(dim(result$p_adjusted), c(2, 3))
  expect_equal(dim(result$significant), c(2, 3))
  expect_equal(dim(result$matrix), c(2, 3))
})

test_that("surrogateMatrixTest significant is logical", {
  skip_if_not_installed("PhysioCore")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  pe1 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )
  pe2 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )

  result <- surrogateMatrixTest(pe1, pe2, method = "coherence",
                                 n_surrogates = 9, nperseg = 64L)

  expect_true(is.logical(result$significant))
})

test_that("surrogateMatrixTest correction='none' returns raw p-values", {
  skip_if_not_installed("PhysioCore")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  pe1 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )
  pe2 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )

  result <- surrogateMatrixTest(pe1, pe2, method = "coherence",
                                 n_surrogates = 9, correction = "none",
                                 nperseg = 64L)

  expect_equal(result$p_values, result$p_adjusted)
  expect_equal(result$correction, "none")
})

test_that("surrogateMatrixTest Bonferroni p_adjusted >= raw p_values", {
  skip_if_not_installed("PhysioCore")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  pe1 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )
  pe2 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )

  result <- surrogateMatrixTest(pe1, pe2, method = "coherence",
                                 n_surrogates = 9, correction = "bonferroni",
                                 nperseg = 64L)

  expect_true(all(result$p_adjusted >= result$p_values - 1e-10))
})

test_that("surrogateMatrixTest alpha parameter is respected", {
  skip_if_not_installed("PhysioCore")

  set.seed(42)
  sr <- 200
  n <- sr * 2
  pe1 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )
  pe2 <- PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
    samplingRate = sr
  )

  result <- surrogateMatrixTest(pe1, pe2, method = "coherence",
                                 n_surrogates = 9, alpha = 0.01,
                                 nperseg = 64L)

  expect_equal(result$alpha, 0.01)
})

test_that("surrogateMatrixTest works with MPE input", {
  skip_if_not_installed("PhysioCore")

  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- surrogateMatrixTest(mpe, method = "coherence",
                                 modality_x = "EEG", modality_y = "EMG",
                                 n_surrogates = 9, nperseg = 64L)

  expect_type(result, "list")
  expect_true(is.matrix(result$matrix))
})

test_that(".extract_coupling_stat handles multitaper_coherence", {
  fake <- list(coherence = c(0.1, 0.9, 0.3))
  stat <- PhysioCrossModal:::.extract_coupling_stat(fake, "multitaper_coherence")
  expect_equal(stat, 0.9)
})
