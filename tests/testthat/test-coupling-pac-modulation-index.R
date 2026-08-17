library(testthat)
library(PhysioCrossModal)

.ws916_pac_signal <- function(n = 2400L, sr = 500, coupled = TRUE,
                              seed = 21L) {
  set.seed(seed)
  time <- (0:(n - 1L)) / sr
  theta <- suppressWarnings(
    PhysioCrossModal:::.bandpass_filter(
      stats::rnorm(n), 5, 7, sr, clamp = FALSE
    )
  )
  theta <- theta / stats::sd(theta)
  phase <- PhysioCrossModal:::.hilbert_phase(theta)
  modulation <- if (coupled) (1 + cos(phase - pi)) / 2 else {
    stats::runif(n)
  }
  theta + (0.1 + modulation) * sin(2 * pi * 50 * time) +
    stats::rnorm(n, sd = 0.25)
}

test_that("modulationIndex observed grids retain comodulogram compatibility", {
  signal <- .ws916_pac_signal()
  for (method in c("tort", "canolty", "ozkurt")) {
    legacy <- comodulogram(
      signal, sr = 500, phase_freqs = c(6, 10), amp_freqs = c(50, 80),
      method = method, phase_bw = 1, amp_bw = 8
    )
    result <- modulationIndex(
      signal, c(6, 10), c(50, 80), sr = 500, method = method,
      phase_bw = 1, amp_bw = 8, n_surrogates = 0
    )
    expect_equal(result$matrix, legacy$matrix, tolerance = 1e-12)
    expect_equal(result$peak, legacy$peak, tolerance = 1e-12)
    expect_s3_class(result, "pac_comodulogram")
    expect_null(result$p_value)
    expect_match(result$warnings, "inference_disabled", all = FALSE)
  }
})

test_that("modulationIndex records the declared grid and realized bands", {
  result <- modulationIndex(
    .ws916_pac_signal(), c(6, 10), c(50, 80), sr = 500,
    phase_bw = 1, amp_bw = 8, n_surrogates = 0
  )
  expect_equal(dim(result$matrix), c(2L, 2L))
  expect_identical(dimnames(result$matrix), list(c("6", "10"), c("50", "80")))
  expect_equal(
    unname(result$realized_bands$phase),
    rbind(c(6, 5, 7, 5, 7), c(10, 9, 11, 9, 11))
  )
  expect_equal(
    unname(result$realized_bands$amplitude),
    rbind(c(50, 42, 58, 42, 58), c(80, 72, 88, 72, 88))
  )
  expect_equal(result$peak$phase_freq,
               result$phase_freqs[which.max(apply(result$matrix, 1, max))])
  expect_true(all(is.finite(result$matrix)))
  expect_true(all(result$matrix >= -sqrt(.Machine$double.eps)))
})

test_that("numerical boundary adjustment is explicit and legacy-compatible", {
  signal <- .ws916_pac_signal()
  legacy <- comodulogram(
    signal, sr = 500, phase_freqs = c(4, 6), amp_freqs = c(50, 80)
  )
  result <- modulationIndex(
    signal, c(4, 6), c(50, 80), sr = 500, n_surrogates = 0
  )
  expect_equal(result$matrix, legacy$matrix, tolerance = 1e-12)
  expect_equal(result$realized_bands$phase["4", "requested_low"], 2)
  expect_equal(result$realized_bands$phase["4", "low"], 2.5)
  expect_match(result$warnings, "requested_band_adjusted", all = FALSE)
})

test_that("PAC formulas match direct definitions", {
  set.seed(916)
  phase <- stats::runif(997, -pi, pi)
  amplitude <- stats::rexp(997) + 0.1
  bins <- 18L

  breaks <- seq(-pi, pi, length.out = bins + 1L)
  index <- cut(phase, breaks, include.lowest = TRUE, labels = FALSE)
  means <- vapply(seq_len(bins), function(i) {
    values <- amplitude[index == i]
    if (length(values)) mean(values) else 0
  }, numeric(1))
  probability <- means / sum(means)
  nonzero <- probability > 0
  tort_reference <- (
    log(bins) + sum(probability[nonzero] * log(probability[nonzero]))
  ) / log(bins)
  canolty_reference <- Mod(sum(amplitude * exp(1i * phase)) /
                             length(amplitude))
  ozkurt_reference <- Mod(sum(amplitude * exp(1i * phase))) /
    sqrt(length(amplitude) * sum(amplitude^2))

  expect_equal(PhysioCrossModal:::.pac_tort(phase, amplitude, bins)$mi,
               tort_reference, tolerance = 1e-12)
  expect_equal(PhysioCrossModal:::.pac_canolty(phase, amplitude),
               canolty_reference, tolerance = 1e-12)
  expect_equal(PhysioCrossModal:::.pac_ozkurt(phase, amplitude),
               ozkurt_reference, tolerance = 1e-12)

  sparse_phase <- c(-3, -2.9, 2.9, 3)
  sparse_amplitude <- c(1, 2, 3, 4)
  expect_equal(
    PhysioCrossModal:::.pac_tort(sparse_phase, sparse_amplitude, bins)$mi,
    local({
      sparse_index <- cut(sparse_phase, breaks, include.lowest = TRUE,
                          labels = FALSE)
      sparse_mean <- vapply(seq_len(bins), function(i) {
        values <- sparse_amplitude[sparse_index == i]
        if (length(values)) mean(values) else 0
      }, numeric(1))
      sparse_probability <- sparse_mean / sum(sparse_mean)
      occupied <- sparse_probability > 0
      (log(bins) + sum(sparse_probability[occupied] *
                         log(sparse_probability[occupied]))) / log(bins)
    }),
    tolerance = 1e-15
  )
})

test_that("pinned Tensorpac formula fixture is intact and matches offline", {
  fixture <- system.file(
    "validation", "fixtures", "tensorpac-0.6.5",
    package = "PhysioCrossModal"
  )
  expect_true(nzchar(fixture))
  checksum <- strsplit(
    readLines(file.path(fixture, "SHA256SUMS")), "  ", fixed = TRUE
  )
  expect_true(all(vapply(checksum, function(parts) {
    path <- normalizePath(file.path(fixture, parts[[2L]]), mustWork = TRUE)
    identical(
      digest::digest(file = path, algo = "sha256", serialize = FALSE),
      parts[[1L]]
    )
  }, logical(1))))

  input <- utils::read.csv(file.path(fixture, "input.csv"))
  expected <- utils::read.csv(file.path(fixture, "expected.csv"))
  observed <- c(
    tort = PhysioCrossModal:::.pac_tort(input$phase, input$amplitude, 18L)$mi,
    canolty = PhysioCrossModal:::.pac_canolty(input$phase, input$amplitude),
    ozkurt = PhysioCrossModal:::.pac_ozkurt(input$phase, input$amplitude)
  )
  reference <- setNames(expected$value, expected$metric)
  expect_equal(observed[names(reference)], reference, tolerance = 1e-12)
})

test_that("coupled jittered PAC peaks and survives grid-wide BH adjustment", {
  signal <- .ws916_pac_signal(n = 4000L, coupled = TRUE)
  old <- options(PhysioCrossModal.pac.keep_surrogates = TRUE)
  on.exit(options(old), add = TRUE)
  result <- modulationIndex(
    signal, c(6, 10), c(50, 80), sr = 500,
    phase_bw = 1, amp_bw = 8, n_surrogates = 99,
    surrogate_type = "timeshift", seed = 301
  )

  expect_equal(result$peak$phase_freq, 6)
  expect_equal(result$peak$amp_freq, 50)
  expect_true(result$significant["6", "50"])
  expect_equal(min(result$p_value), 1 / 100)
  expect_true(all(result$p_value > 0))

  null <- result$surrogate_array
  null_flat <- matrix(null, nrow = length(result$matrix))
  manual_p <- matrix(
    (rowSums(null_flat >= as.vector(result$matrix)) + 1) / 100,
    nrow = 2, dimnames = dimnames(result$matrix)
  )
  manual_bh <- matrix(
    stats::p.adjust(as.vector(manual_p), method = "BH"),
    nrow = 2, dimnames = dimnames(result$matrix)
  )
  expect_equal(result$p_value, manual_p)
  expect_equal(result$p_adjusted, manual_bh)
  transposed <- t(unname(manual_bh))
  dimnames(transposed) <- dimnames(manual_bh)
  expect_false(isTRUE(all.equal(result$p_adjusted, transposed)))

  set.seed(301)
  first_surrogate <- PhysioCrossModal:::.timeshift_surrogate(signal)
  phase_bands <- result$realized_bands$phase
  amp_bands <- result$realized_bands$amplitude
  phases <- PhysioCrossModal:::.pac_filter_phases(
    signal, phase_bands, 500, 4L
  )
  amplitudes <- PhysioCrossModal:::.pac_filter_amplitudes(
    first_surrogate, amp_bands, 500, 4L
  )
  names(amplitudes) <- rownames(amp_bands)
  first_reference <- PhysioCrossModal:::.pac_grid_from_components(
    phases, amplitudes, "tort", 18L, 500, phase_bands, 4L
  )
  expect_equal(null[, , 1], first_reference, tolerance = 1e-12)
})

test_that("seeded sequential and parallel inference are identical and restore RNG", {
  signal <- .ws916_pac_signal(n = 1600L)
  set.seed(700)
  before <- .Random.seed
  sequential <- modulationIndex(
    signal, 6, 50, sr = 500, phase_bw = 1, amp_bw = 8,
    n_surrogates = 20, seed = 51, cores = 1
  )
  expect_identical(.Random.seed, before)
  parallel <- modulationIndex(
    signal, 6, 50, sr = 500, phase_bw = 1, amp_bw = 8,
    n_surrogates = 20, seed = 51, cores = 2
  )
  expect_identical(.Random.seed, before)
  expect_equal(parallel$matrix, sequential$matrix)
  expect_equal(parallel$p_value, sequential$p_value)
  expect_equal(parallel$threshold, sequential$threshold)

  disabled <- modulationIndex(
    signal, 6, 50, sr = 500, phase_bw = 1, amp_bw = 8,
    n_surrogates = 0, seed = 52
  )
  expect_identical(.Random.seed, before)
  expect_null(disabled$threshold)
})

test_that("modulationIndex restores an initially absent RNG state", {
  signal <- .ws916_pac_signal(n = 1200L)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", saved_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  invisible(modulationIndex(
    signal, 6, 50, sr = 500, phase_bw = 1, amp_bw = 8,
    n_surrogates = 20, seed = 52
  ))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("modulationIndex fails before partial computation on invalid inputs", {
  signal <- .ws916_pac_signal(n = 500L)
  base <- list(x = signal, phase_freqs = c(6, 10), amp_freqs = c(50, 80),
               sr = 500, phase_bw = 1, amp_bw = 8, n_surrogates = 0)
  call_mi <- function(...) do.call(modulationIndex, modifyList(base, list(...)))

  expect_error(call_mi(phase_freqs = c(10, 6)), "strictly increasing")
  expect_error(call_mi(phase_freqs = c(6, 6)), "unique")
  expect_error(call_mi(phase_freqs = c(0.5, 6)), "strictly within")
  expect_error(call_mi(amp_freqs = c(50, 245)), "strictly within")
  expect_error(call_mi(phase_bw = 0), "phase_bw")
  expect_error(call_mi(amp_bw = Inf), "amplitude_bw")
  expect_error(call_mi(n_bins = 1), "n_bins")
  expect_error(call_mi(order = 2.5), "order")
  expect_error(call_mi(n_surrogates = 19), "at least 20")
  expect_error(call_mi(n_surrogates = 20, p_adjust_method = "invalid"),
               "p_adjust_method")
  expect_error(call_mi(n_surrogates = 20, seed = 1.5), "seed")
  expect_error(call_mi(n_surrogates = 20, cores = 0), "cores")
  expect_error(call_mi(channels_x = c(1, 2)), "channels_x")
  expect_error(call_mi(x = matrix(signal, ncol = 1)), "must be a vector")
  expect_error(call_mi(x = rep(1, 500)), "non-constant")
  expect_error(call_mi(x = c(signal[-1], NA_real_)), "finite")
  expect_error(call_mi(x = signal[1:20], phase_freqs = 6, amp_freqs = 50),
               "too short")
  expect_error(call_mi(y = signal[-1]), "equal lengths")
})

test_that("plotComodulogram validates and applies significance masks", {
  skip_if_not_installed("ggplot2")
  signal <- .ws916_pac_signal(n = 1200L)
  legacy <- comodulogram(
    signal, sr = 500, phase_freqs = c(6, 10),
    amp_freqs = c(50, 70, 90), phase_bw = 1, amp_bw = 8
  )
  expect_s3_class(plotComodulogram(legacy), "ggplot")

  result <- modulationIndex(
    signal, c(6, 10), c(50, 70, 90), sr = 500,
    phase_bw = 1, amp_bw = 8, n_surrogates = 20, seed = 13
  )
  automatic <- plotComodulogram(result, nonsignificant_alpha = 0.15)
  built <- ggplot2::ggplot_build(automatic)
  expect_equal(
    sort(unique(built$data[[1]]$alpha)),
    sort(unique(ifelse(as.vector(result$significant), 1, 0.15)))
  )

  logical_mask <- result$significant
  numeric_mask <- result$p_adjusted
  expect_s3_class(plotComodulogram(result, mask = logical_mask), "ggplot")
  expect_s3_class(plotComodulogram(result, mask = numeric_mask), "ggplot")
  expect_error(plotComodulogram(result, mask = t(logical_mask)),
               "dimensions and dimnames")
  bad_names <- logical_mask
  dimnames(bad_names)[[1]] <- rev(dimnames(bad_names)[[1]])
  expect_error(plotComodulogram(result, mask = bad_names), "dimnames")
  logical_mask[1] <- NA
  expect_error(plotComodulogram(result, mask = logical_mask), "cannot contain")
  numeric_mask[1] <- Inf
  expect_error(plotComodulogram(result, mask = numeric_mask), "finite p-values")
  expect_error(plotComodulogram(result, mask = TRUE), "matrix")
})
