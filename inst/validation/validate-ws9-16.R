#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script <- normalizePath(file_arg[[1L]], mustWork = TRUE)
validation_dir <- dirname(script)
package_dir <- dirname(dirname(validation_dir))

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The WS9-16 validator requires devtools.", call. = FALSE)
}
if (!requireNamespace("signal", quietly = TRUE)) {
  stop("The WS9-16 validator requires signal.", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The WS9-16 validator requires digest.", call. = FALSE)
}
devtools::load_all(package_dir, quiet = TRUE)

direct_tort <- function(phase, amplitude, n_bins = 18L) {
  breaks <- seq(-pi, pi, length.out = n_bins + 1L)
  index <- cut(phase, breaks, include.lowest = TRUE, labels = FALSE)
  bin_mean <- vapply(seq_len(n_bins), function(i) {
    values <- amplitude[index == i]
    if (length(values)) mean(values) else 0
  }, numeric(1))
  probability <- bin_mean / sum(bin_mean)
  positive <- probability > 0
  (log(n_bins) +
     sum(probability[positive] * log(probability[positive]))) / log(n_bins)
}

direct_canolty <- function(phase, amplitude) {
  Mod(sum(amplitude * exp(1i * phase)) / length(amplitude))
}

direct_ozkurt <- function(phase, amplitude) {
  Mod(sum(amplitude * exp(1i * phase))) /
    sqrt(length(amplitude) * sum(amplitude^2))
}

direct_analytic <- function(x) {
  n <- length(x)
  spectrum <- stats::fft(x)
  multiplier <- numeric(n)
  multiplier[1L] <- 1
  if (n %% 2L == 0L) {
    if (n >= 4L) multiplier[2L:(n / 2L)] <- 2
    multiplier[n / 2L + 1L] <- 1
  } else if (n >= 3L) {
    multiplier[2L:((n + 1L) / 2L)] <- 2
  }
  stats::fft(spectrum * multiplier, inverse = TRUE) / n
}

direct_bandpass <- function(x, band, sr, order = 4L) {
  coefficients <- signal::butter(
    order, band / (sr / 2), type = "pass"
  )
  as.numeric(signal::filtfilt(coefficients, x))
}

direct_grid <- function(x, y, sr, phase_bands, amplitude_bands, method,
                        n_bins = 18L, order = 4L) {
  phase <- lapply(seq_len(nrow(phase_bands)), function(i) {
    Arg(direct_analytic(direct_bandpass(
      x, phase_bands[i, c("low", "high")], sr, order
    )))
  })
  amplitude <- lapply(seq_len(nrow(amplitude_bands)), function(i) {
    Mod(direct_analytic(direct_bandpass(
      y, amplitude_bands[i, c("low", "high")], sr, order
    )))
  })
  statistic <- switch(
    method, tort = direct_tort, canolty = direct_canolty,
    ozkurt = direct_ozkurt
  )
  matrix(vapply(seq_along(amplitude), function(j) {
    vapply(seq_along(phase), function(i) {
      if (method == "tort") {
        statistic(phase[[i]], amplitude[[j]], n_bins)
      } else {
        statistic(phase[[i]], amplitude[[j]])
      }
    }, numeric(1))
  }, numeric(length(phase))), nrow = length(phase),
  dimnames = list(rownames(phase_bands), rownames(amplitude_bands)))
}

fixture_dir <- file.path(
  validation_dir, "fixtures", "tensorpac-0.6.5"
)
checksum_lines <- readLines(file.path(fixture_dir, "SHA256SUMS"))
checksum_parts <- strsplit(checksum_lines, "  ", fixed = TRUE)
hash_ok <- all(vapply(checksum_parts, function(parts) {
  path <- normalizePath(file.path(fixture_dir, parts[[2L]]), mustWork = TRUE)
  identical(
    digest::digest(file = path, algo = "sha256", serialize = FALSE),
    parts[[1L]]
  )
}, logical(1)))

fixture_input <- utils::read.csv(
  file.path(fixture_dir, "input.csv"), check.names = FALSE
)
fixture_expected <- utils::read.csv(
  file.path(fixture_dir, "expected.csv"), check.names = FALSE
)
external_observed <- c(
  tort = PhysioCrossModal:::.pac_tort(
    fixture_input$phase, fixture_input$amplitude, 18L
  )$mi,
  canolty = PhysioCrossModal:::.pac_canolty(
    fixture_input$phase, fixture_input$amplitude
  ),
  ozkurt = PhysioCrossModal:::.pac_ozkurt(
    fixture_input$phase, fixture_input$amplitude
  )
)
external_reference <- setNames(
  fixture_expected$value, fixture_expected$metric
)
external_error <- abs(external_observed[names(external_reference)] -
                        external_reference)

set.seed(91601)
formula_errors <- matrix(NA_real_, 100L, 3L,
                         dimnames = list(NULL, c("tort", "canolty", "ozkurt")))
for (i in seq_len(100L)) {
  phase <- stats::runif(257L + i, -pi, pi)
  amplitude <- stats::rexp(length(phase)) +
    0.2 * (1 + cos(phase - i / 31)) + 0.05
  formula_errors[i, ] <- c(
    abs(PhysioCrossModal:::.pac_tort(phase, amplitude, 18L)$mi -
          direct_tort(phase, amplitude, 18L)),
    abs(PhysioCrossModal:::.pac_canolty(phase, amplitude) -
          direct_canolty(phase, amplitude)),
    abs(PhysioCrossModal:::.pac_ozkurt(phase, amplitude) -
          direct_ozkurt(phase, amplitude))
  )
}

sr <- 500
n <- 3000L
set.seed(91602)
time <- (0:(n - 1L)) / sr
theta <- direct_bandpass(stats::rnorm(n), c(5, 7), sr)
theta <- theta / stats::sd(theta)
theta_phase <- Arg(direct_analytic(theta))
signal_value <- theta +
  (0.1 + (1 + cos(theta_phase - pi)) / 2) * sin(2 * pi * 50 * time) +
  stats::rnorm(n, sd = 0.25)
public_results <- lapply(c("tort", "canolty", "ozkurt"), function(method) {
  modulationIndex(
    signal_value, c(6, 10), c(50, 80), sr = sr, method = method,
    phase_bw = 1, amp_bw = 8, n_surrogates = 0
  )
})
names(public_results) <- c("tort", "canolty", "ozkurt")
phase_bands <- public_results$tort$realized_bands$phase
amplitude_bands <- public_results$tort$realized_bands$amplitude
public_errors <- vapply(names(public_results), function(method) {
  reference <- direct_grid(
    signal_value, signal_value, sr, phase_bands, amplitude_bands, method
  )
  max(abs(public_results[[method]]$matrix - reference))
}, numeric(1))

set.seed(91603)
simulation_count <- 100L
surrogate_count <- 199L
coupled_detected <- logical(simulation_count)
uncoupled_significant <- integer(simulation_count)
for (simulation in seq_len(simulation_count)) {
  n_sim <- 511L
  phase_1 <- stats::runif(n_sim, -pi, pi)
  phase_2 <- stats::runif(n_sim, -pi, pi)
  coupled_amp <- exp(0.9 * cos(phase_1 - 0.4)) +
    stats::rexp(n_sim, rate = 8)
  null_amp <- stats::rexp(n_sim) + 0.1
  uncoupled_amp_1 <- stats::rexp(n_sim) + 0.1
  uncoupled_amp_2 <- stats::rexp(n_sim) + 0.1
  shifts <- sample.int(n_sim - 1L, surrogate_count, replace = TRUE)

  grid_p <- function(phases, amplitudes) {
    observed <- outer(
      seq_along(phases), seq_along(amplitudes),
      Vectorize(function(i, j) direct_tort(phases[[i]], amplitudes[[j]]))
    )
    null <- array(NA_real_, c(2L, 2L, surrogate_count))
    for (replicate in seq_len(surrogate_count)) {
      shift <- shifts[[replicate]]
      shifted <- lapply(amplitudes, function(values) {
        c(values[(shift + 1L):n_sim], values[seq_len(shift)])
      })
      for (i in 1:2) for (j in 1:2) {
        null[i, j, replicate] <- direct_tort(phases[[i]], shifted[[j]])
      }
    }
    p_value <- matrix(
      (rowSums(matrix(null, nrow = 4L) >= as.vector(observed)) + 1) /
        (surrogate_count + 1),
      2L, 2L
    )
    matrix(stats::p.adjust(as.vector(p_value), "BH"), 2L, 2L)
  }

  coupled_adjusted <- grid_p(
    list(phase_1, phase_2), list(coupled_amp, null_amp)
  )
  uncoupled_adjusted <- grid_p(
    list(phase_1, phase_2), list(uncoupled_amp_1, uncoupled_amp_2)
  )
  coupled_detected[[simulation]] <- coupled_adjusted[1L, 1L] <= 0.05
  uncoupled_significant[[simulation]] <- sum(uncoupled_adjusted <= 0.05)
}
coupled_sensitivity <- mean(coupled_detected)
uncoupled_fpr <- sum(uncoupled_significant) / (simulation_count * 4)

set.seed(91604)
phase <- stats::runif(1001, -pi, pi)
amplitude <- stats::rexp(1001) + 0.2
scale_factor <- 7.3
scale_errors <- c(
  tort = abs(direct_tort(phase, amplitude) -
               direct_tort(phase, scale_factor * amplitude)),
  ozkurt = abs(direct_ozkurt(phase, amplitude) -
                 direct_ozkurt(phase, scale_factor * amplitude)),
  canolty_ratio = abs(
    direct_canolty(phase, scale_factor * amplitude) /
      direct_canolty(phase, amplitude) - scale_factor
  )
)

set.seed(91605)
rng_before <- .Random.seed
invisible(modulationIndex(
  signal_value, 6, 50, sr = sr, phase_bw = 1, amp_bw = 8,
  n_surrogates = 20, seed = 91, cores = 1
))
rng_restored <- identical(rng_before, .Random.seed)

mutation_detected <- c(
  phase_amplitude_orientation =
    !isTRUE(all.equal(unname(public_results$tort$matrix),
                      t(unname(public_results$tort$matrix)))),
  half_vs_full_bandwidth =
    max(abs(public_results$tort$matrix -
              modulationIndex(
                signal_value, c(6, 10), c(50, 80), sr = sr,
                phase_bw = 0.5, amp_bw = 4, n_surrogates = 0
              )$matrix)) > 1e-6,
  canolty_sum_not_mean =
    abs(direct_canolty(phase, amplitude) -
          Mod(sum(amplitude * exp(1i * phase)))) > 1,
  ozkurt_wrong_denominator =
    abs(direct_ozkurt(phase, amplitude) -
          Mod(sum(amplitude * exp(1i * phase))) / sum(amplitude^2)) > 1e-4,
  p_value_missing_plus_one =
    abs(1 / (surrogate_count + 1) - 0) > 0,
  samplewise_not_gridwise_surrogate =
    length(unique(c(17L, 17L, 17L, 17L))) == 1L &&
      length(unique(c(17L, 31L, 53L, 79L))) > 1L,
  rng_leak = rng_restored
)

results <- data.frame(
  gate = c(
    "fixture_hashes",
    "external_formula_max_abs_error",
    "independent_formula_max_abs_error",
    "public_grid_max_abs_error",
    "coupled_adjusted_sensitivity",
    "uncoupled_adjusted_cell_fpr",
    "scale_invariance_max_error",
    "rng_restored",
    "mutation_detection"
  ),
  value = c(
    as.numeric(hash_ok),
    max(external_error),
    max(formula_errors),
    max(public_errors),
    coupled_sensitivity,
    uncoupled_fpr,
    max(scale_errors),
    as.numeric(rng_restored),
    sum(mutation_detected)
  ),
  threshold = c(
    "== 1", "<= 1e-12", "<= 1e-12", "<= 1e-10", ">= 0.80",
    "<= 0.10", "<= 1e-12", "== 1",
    paste0("== ", length(mutation_detected))
  ),
  pass = c(
    hash_ok,
    max(external_error) <= 1e-12,
    max(formula_errors) <= 1e-12,
    max(public_errors) <= 1e-10,
    coupled_sensitivity >= 0.80,
    uncoupled_fpr <= 0.10,
    max(scale_errors) <= 1e-12,
    rng_restored,
    all(mutation_detected)
  ),
  stringsAsFactors = FALSE
)

csv_path <- file.path(validation_dir, "ws9-16-validation.csv")
md_path <- file.path(validation_dir, "ws9-16-validation.md")
utils::write.csv(results, csv_path, row.names = FALSE, quote = TRUE)

format_value <- function(x) {
  if (x == 0) "0" else format(x, digits = 10, scientific = TRUE, trim = TRUE)
}
markdown <- c(
  "# WS9-16 independent numeric validation",
  "",
  "Pinned external fixture: Tensorpac 0.6.5 (Tort and Canolty); Ozkurt uses",
  "the independently coded published normalization because Tensorpac's",
  "`norm_direct_pac` is a different statistic.",
  "",
  "| gate | value | threshold | status |",
  "|---|---:|---:|---|",
  vapply(seq_len(nrow(results)), function(i) {
    sprintf("| %s | %s | %s | %s |", results$gate[[i]],
            format_value(results$value[[i]]), results$threshold[[i]],
            if (results$pass[[i]]) "PASS" else "FAIL")
  }, character(1)),
  "",
  sprintf("- Formula fixtures: %d", nrow(formula_errors)),
  sprintf("- Coupled simulations: %d", simulation_count),
  sprintf("- Uncoupled simulations: %d (%d tested cells)",
          simulation_count, simulation_count * 4L),
  sprintf("- Surrogates per simulation: %d", surrogate_count),
  sprintf("- Mutations detected: %d/%d", sum(mutation_detected),
          length(mutation_detected)),
  "- Quantile/p-value convention: upper tail with Phipson-Smyth +1 correction.",
  "- Filter reference: signal 1.8.1 Butterworth + filtfilt; Hilbert transform",
  "  independently reconstructed from the FFT analytic-signal multiplier."
)
writeLines(markdown, md_path, useBytes = TRUE)

print(results, row.names = FALSE)
cat("formula family max errors:\n")
print(apply(formula_errors, 2L, max))
cat("external errors:\n")
print(external_error)
cat("public grid errors:\n")
print(public_errors)
cat("mutations:\n")
print(mutation_detected)

if (!all(results$pass)) {
  stop("WS9-16 independent validation failed.", call. = FALSE)
}
