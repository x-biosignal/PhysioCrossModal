# PhysioCrossModal now takes the ecosystem's one multimodal container as input
# instead of defining a rival class. These tests cover that connection and the
# two defects it was blocking.

mk <- function(n, sr, seed = 1) {
  set.seed(seed)
  PhysioExperiment(assays = list(raw = matrix(rnorm(n * 2), n, 2)),
                   samplingRate = sr)
}

test_that("this package defines no container class of its own", {
  expect_equal(length(methods::getClasses(asNamespace("PhysioCrossModal"))), 0L)
  m <- MultiPhysioExperiment(a = mk(50, 100))
  expect_identical(attr(class(m), "package"), "PhysioExperiment")
})

test_that("the legacy accessors keep their names and return types", {
  m <- MultiPhysioExperiment(EEG = mk(500, 250), EMG = mk(500, 250))
  expect_type(experiments(m), "list")
  expect_equal(modalities(m), c("EEG", "EMG"))
  expect_equal(unname(samplingRates(m)), c(250, 250))
  expect_equal(nModalities(m), 2L)
  expect_equal(length(m), 2L)
  expect_equal(names(m), c("EEG", "EMG"))
  expect_s4_class(alignment(m), "DataFrame")
})

# ---- defect 1: pair extraction ignored the clock ----------------------------

test_that("simultaneous streams give exactly the result they always did", {
  # Equal offsets must take the historical path untouched, so existing numbers
  # on simultaneously started data are preserved.
  m <- MultiPhysioExperiment(streams = list(EEG = mk(512, 256, 1),
                                            EMG = mk(512, 256, 2)),
                             offsets = c(EEG = 0, EMG = 0))
  pair <- PhysioCrossModal:::.extract_signal_pair(
    m, NULL, NULL, modality_x = "EEG", modality_y = "EMG",
    channels_x = 1L, channels_y = 1L)

  expect_true(pair$alignment$simultaneous_start)
  expect_equal(pair$x, as.numeric(SummarizedExperiment::assay(
    experiments(m)$EEG, 1L)[, 1]))
  expect_equal(pair$y, as.numeric(SummarizedExperiment::assay(
    experiments(m)$EMG, 1L)[, 1]))
})

test_that("streams that start at different times are compared on their overlap", {
  # EMG starts 1 s late. Head-to-head comparison would pair EEG at t=0 with EMG
  # at t=1 -- a one-second error that produces a plausible-looking number.
  # 400 samples at 100 Hz spans (400-1)/100 = 3.99 s, not 4 s
  eeg <- mk(400, 100, 1)   # t = 0.00 .. 3.99
  emg <- mk(400, 100, 2)   # t = 1.00 .. 4.99
  m <- MultiPhysioExperiment(streams = list(EEG = eeg, EMG = emg),
                             offsets = c(EEG = 0, EMG = 1))
  pair <- PhysioCrossModal:::.extract_signal_pair(
    m, NULL, NULL, modality_x = "EEG", modality_y = "EMG",
    channels_x = 1L, channels_y = 1L)

  expect_false(pair$alignment$simultaneous_start)
  expect_equal(pair$alignment$overlap, c(1, 3.99))
  # the overlap is 2.99 s at 100 Hz -> 300 samples of each
  expect_equal(length(pair$x), 300L)
  expect_equal(length(pair$y), 300L)
  # x starts at its sample for t = 1 s, y at its own first sample (t = 1 s)
  eeg_raw <- as.numeric(SummarizedExperiment::assay(eeg, 1L)[, 1])
  emg_raw <- as.numeric(SummarizedExperiment::assay(emg, 1L)[, 1])
  expect_equal(pair$x[1], eeg_raw[101])
  expect_equal(pair$y[1], emg_raw[1])
  # and the pairing is time-correct all the way along
  expect_equal(pair$x, eeg_raw[101:400])
})

test_that("streams that never overlap fail loudly instead of returning a number", {
  m <- MultiPhysioExperiment(streams = list(EEG = mk(100, 100, 1),
                                            EMG = mk(100, 100, 2)),
                             offsets = c(EEG = 0, EMG = 50))
  expect_error(
    PhysioCrossModal:::.extract_signal_pair(
      m, NULL, NULL, modality_x = "EEG", modality_y = "EMG",
      channels_x = 1L, channels_y = 1L),
    "do not overlap in time")
})

test_that("the alignment actually applied is reported with the pair", {
  m <- MultiPhysioExperiment(streams = list(EEG = mk(200, 100, 1),
                                            EMG = mk(400, 200, 2)),
                             offsets = c(EEG = 0, EMG = 0.5))
  a <- PhysioCrossModal:::.extract_signal_pair(
    m, NULL, NULL, modality_x = "EEG", modality_y = "EMG",
    channels_x = 1L, channels_y = 1L)$alignment
  expect_equal(a$output_rate, 100)
  expect_true(a$resampled)
  expect_match(a$method, "interpolated \\(linear\\)")
  expect_equal(a$offsets, c(0, 0.5))
  expect_gt(a$n_samples, 0L)
})

# ---- defect 2: the cache key ignored the signals and the clock --------------

test_that("the cache is reused for identical inputs", {
  m <- MultiPhysioExperiment(EEG = mk(512, 256, 1), EMG = mk(512, 256, 2))
  couplingResults(m) <- NULL
  a <- couplingAnalysisCached(m, "coherence", "EEG", "EMG")
  b <- couplingAnalysisCached(a$mpe, "coherence", "EEG", "EMG")
  expect_false(a$cached)
  expect_true(b$cached)
})

test_that("changing the SIGNALS invalidates the cache", {
  m <- MultiPhysioExperiment(EEG = mk(512, 256, 1), EMG = mk(512, 256, 2))
  couplingResults(m) <- NULL
  first <- couplingAnalysisCached(m, "coherence", "EEG", "EMG")

  m2 <- MultiPhysioExperiment(EEG = mk(512, 256, 99), EMG = mk(512, 256, 2))
  second <- couplingAnalysisCached(m2, "coherence", "EEG", "EMG")
  expect_false(second$cached)   # the old key would have reported TRUE here
})

test_that("changing the CLOCK invalidates the cache", {
  eeg <- mk(512, 256, 1); emg <- mk(512, 256, 2)
  m <- MultiPhysioExperiment(streams = list(EEG = eeg, EMG = emg),
                             offsets = c(EEG = 0, EMG = 0))
  couplingResults(m) <- NULL
  invisible(couplingAnalysisCached(m, "coherence", "EEG", "EMG"))

  shifted <- MultiPhysioExperiment(streams = list(EEG = eeg, EMG = emg),
                                   offsets = c(EEG = 0, EMG = 0.25))
  out <- couplingAnalysisCached(shifted, "coherence", "EEG", "EMG")
  expect_false(out$cached)
})

test_that("changing the method or its parameters invalidates the cache", {
  m <- MultiPhysioExperiment(EEG = mk(512, 256, 1), EMG = mk(512, 256, 2))
  couplingResults(m) <- NULL
  invisible(couplingAnalysisCached(m, "coherence", "EEG", "EMG"))
  expect_false(couplingAnalysisCached(m, "crosscorrelation", "EEG", "EMG")$cached)
})

test_that("assigned results are bound to the container they were assigned for", {
  m <- MultiPhysioExperiment(EEG = mk(256, 128, 1), EMG = mk(256, 128, 2))
  couplingResults(m) <- list(a = 1)
  expect_identical(couplingResults(m)$a, 1)

  # the same assignment does not leak onto a container with different signals
  other <- MultiPhysioExperiment(EEG = mk(256, 128, 77), EMG = mk(256, 128, 2))
  expect_equal(length(couplingResults(other)), 0L)

  couplingResults(m) <- NULL
  expect_equal(length(couplingResults(m)), 0L)
})

# ---- regression for the 2026-09-26 review ----------------------------------

test_that("a start difference smaller than one sample is removed, not carried", {
  # each sample's value IS its own true time, so a correct alignment pairs
  # equal values. The offset is half a sample interval at 100 Hz.
  tmk <- function(times) PhysioExperiment(
    assays = list(raw = matrix(times, ncol = 1)), samplingRate = 100)
  a <- tmk(seq(0, 0.99, by = 0.01))
  b <- tmk(seq(0.005, 0.995, by = 0.01))
  m <- MultiPhysioExperiment(streams = list(a = a, b = b),
                             offsets = c(a = 0, b = 0.005))
  p <- PhysioCrossModal:::.extract_signal_pair(
    m, NULL, NULL, modality_x = "a", modality_y = "b",
    channels_x = 1L, channels_y = 1L)

  expect_lt(max(abs(p$x - p$y)), 1e-9)   # head-to-head pairing left 0.005
  expect_false(p$alignment$simultaneous_start)
  expect_match(p$alignment$method, "onto one")
})

test_that("out-of-phase grids at the same rate are aligned exactly", {
  # Same rate means no anti-aliasing filter runs, so the "value is its own
  # time" trick isolates the time alignment itself.
  tmk <- function(times) PhysioExperiment(
    assays = list(raw = matrix(times, ncol = 1)), samplingRate = 100)
  a <- tmk(seq(0, 1, by = 0.01))
  b <- tmk(seq(0.0033, 1.0033, by = 0.01))
  p <- PhysioCrossModal:::.extract_signal_pair(
    MultiPhysioExperiment(streams = list(a = a, b = b),
                          offsets = c(a = 0, b = 0.0033)),
    NULL, NULL, modality_x = "a", modality_y = "b",
    channels_x = 1L, channels_y = 1L)
  expect_lt(max(abs(p$x - p$y)), 1e-9)
  expect_equal(length(p$x), length(p$y))
})

test_that("different rates are placed on one grid and the conversion is reported", {
  # Here an anti-aliasing filter legitimately changes y's values before
  # interpolation, so what is asserted is the grid and the record of it, not
  # value equality.
  tmk <- function(n, sr, off) PhysioExperiment(
    assays = list(raw = matrix(sin(2 * pi * 3 * (off + (seq_len(n) - 1) / sr)), ncol = 1)),
    samplingRate = sr)
  m <- MultiPhysioExperiment(streams = list(a = tmk(101, 100, 0), b = tmk(201, 200, 0.0033)),
                             offsets = c(a = 0, b = 0.0033))
  p <- PhysioCrossModal:::.extract_signal_pair(
    m, NULL, NULL, modality_x = "a", modality_y = "b",
    channels_x = 1L, channels_y = 1L)
  expect_equal(length(p$x), length(p$y))
  expect_equal(p$alignment$output_rate, 100)
  expect_true(p$alignment$resampled)
  expect_equal(unname(p$alignment$grid[["first"]]), 0.0033)
  expect_match(p$alignment$method, "low-pass")
  # a 3 Hz sine survives a 45 Hz anti-alias filter well enough to stay aligned
  expect_lt(max(abs(p$x - p$y)), 0.05)
})

# ---- regressions for the 2026-09-27 assessment ------------------------------
# Both defects survived the previous fix pass because measured sample times were
# made authoritative in the foundation but not consulted here.

tmk_t <- function(times) PhysioExperiment(
  assays = list(raw = matrix(times, ncol = 1)), samplingRate = 100,
  rowData = S4Vectors::DataFrame(time_from_t0 = times))

test_that("an equal start offset is not treated as an equal time correspondence", {
  # both streams begin at t = 0, but `a` has a one-second hole. Pairing element
  # k with element k compared t = 3.57 against t = 2.57.
  a <- tmk_t(c(seq(0, 1.99, by = 0.01), seq(3, 3.99, by = 0.01)))
  b <- tmk_t(seq(0, 2.99, by = 0.01))
  m <- MultiPhysioExperiment(streams = list(a = a, b = b), offsets = c(a = 0, b = 0))

  expect_error(
    PhysioCrossModal:::.extract_signal_pair(m, NULL, NULL, modality_x = "a",
                                            modality_y = "b", channels_x = 1L,
                                            channels_y = 1L),
    "no samples between")
})

test_that("a gap outside the compared window is not an obstacle", {
  a <- tmk_t(c(seq(0, 1.99, by = 0.01), seq(3, 3.99, by = 0.01)))
  b <- tmk_t(seq(0, 1.99, by = 0.01))
  p <- PhysioCrossModal:::.extract_signal_pair(
    MultiPhysioExperiment(streams = list(a = a, b = b), offsets = c(a = 0, b = 0)),
    NULL, NULL, modality_x = "a", modality_y = "b",
    channels_x = 1L, channels_y = 1L)
  expect_lt(max(abs(p$x - p$y)), 1e-9)
  expect_false(p$alignment$regular_common_start)
  expect_true(all(p$alignment$measured_times))
})

test_that("regular streams that start together keep the historical fast path", {
  # no measured times, equal offsets: the two time vectors are provably
  # identical, so the published numbers for this case must not move
  mkr <- function(v) PhysioExperiment(assays = list(raw = matrix(v, ncol = 1)),
                                      samplingRate = 100)
  set.seed(7); a <- mkr(rnorm(300)); set.seed(8); b <- mkr(rnorm(300))
  p <- PhysioCrossModal:::.extract_signal_pair(
    MultiPhysioExperiment(streams = list(a = a, b = b), offsets = c(a = 0, b = 0)),
    NULL, NULL, modality_x = "a", modality_y = "b",
    channels_x = 1L, channels_y = 1L)
  expect_true(p$alignment$regular_common_start)
  expect_equal(p$x, as.numeric(SummarizedExperiment::assay(a, 1L)[, 1]))
  expect_equal(p$y, as.numeric(SummarizedExperiment::assay(b, 1L)[, 1]))
})

test_that("changing only the measured times invalidates the cache", {
  set.seed(1); v1 <- rnorm(400); set.seed(2); v2 <- rnorm(400)
  pe <- function(v, times) PhysioExperiment(
    assays = list(raw = matrix(v, ncol = 1)), samplingRate = 100,
    rowData = S4Vectors::DataFrame(time_from_t0 = times))
  t0 <- seq(0, 3.99, by = 0.01)

  # a half-sample shift of Y: the signal values, rate and offsets are untouched
  m1 <- MultiPhysioExperiment(streams = list(X = pe(v1, t0), Y = pe(v2, t0)),
                              offsets = c(X = 0, Y = 0))
  m2 <- MultiPhysioExperiment(streams = list(X = pe(v1, t0), Y = pe(v2, t0 + 0.005)),
                              offsets = c(X = 0, Y = 0))
  couplingResults(m1) <- NULL; couplingResults(m2) <- NULL

  first <- couplingAnalysisCached(m1, "coherence", "X", "Y")
  after <- couplingAnalysisCached(m2, "coherence", "X", "Y")
  expect_false(after$cached)            # the old key reported TRUE here

  # and what the cache would return equals a fresh computation
  direct <- couplingAnalysis(mpe = m2, method = "coherence",
                             modality_x = "X", modality_y = "Y")
  expect_equal(after$result$coherence, direct$coherence)
  # the stale answer really was a different answer
  expect_gt(max(abs(after$result$coherence - first$result$coherence)), 0.01)
})

test_that("the public result records the inputs it was computed from", {
  mkr <- function(v) PhysioExperiment(assays = list(raw = matrix(v, ncol = 1)),
                                      samplingRate = 100)
  set.seed(3); m <- MultiPhysioExperiment(X = mkr(rnorm(400)), Y = mkr(rnorm(400)))
  out <- couplingAnalysis(mpe = m, method = "coherence",
                          modality_x = "X", modality_y = "Y")
  expect_false(is.null(out$input))
  expect_equal(out$input$modality_x, "X")
  expect_equal(out$input$output_rate, 100)
  expect_false(is.null(out$input$implementation))
})

test_that("the fast-path flag says what it means, not more", {
  # 100 Hz beside 200 Hz with a shared start takes the fast path, and the flag
  # is named for what actually holds: both streams are on reconstructed regular
  # grids from a common origin. Their time vectors are NOT identical.
  mkr <- function(n, sr) PhysioExperiment(
    assays = list(raw = matrix(rnorm(n), ncol = 1)), samplingRate = sr)
  set.seed(1)
  p <- PhysioCrossModal:::.extract_signal_pair(
    MultiPhysioExperiment(streams = list(a = mkr(200, 100), b = mkr(400, 200)),
                          offsets = c(a = 0, b = 0)),
    NULL, NULL, modality_x = "a", modality_y = "b",
    channels_x = 1L, channels_y = 1L)
  expect_true(p$alignment$regular_common_start)
  expect_true(p$alignment$resampled)          # the rates differ, so y was converted
  expect_null(p$alignment$identical_grids)    # the old, overclaiming name is gone
})

test_that("the gap criterion is the adopted 1.5x rule, with stated boundaries", {
  tt <- seq(0, 0.99, by = 0.01)
  # a single dropped sample is 2x the period: flagged
  dropped <- tt[-50]
  expect_error(PhysioCrossModal:::.reject_gap_in_window(dropped, 100, "s", 0, 0.99),
               "no samples between")
  # jitter inside the threshold is not flagged -- the rule separates jitter from
  # a dropout, it does not certify regular sampling
  jittered <- tt + c(0, rep(c(0.002, -0.002), length.out = length(tt) - 1))
  expect_true(PhysioCrossModal:::.reject_gap_in_window(jittered, 100, "s", 0, 0.99))
  # a gap outside the window is not the window's problem
  expect_true(PhysioCrossModal:::.reject_gap_in_window(dropped, 100, "s", 0, 0.3))
  # the factor is a parameter, so the criterion can be stated and varied
  expect_true(PhysioCrossModal:::.reject_gap_in_window(dropped, 100, "s", 0, 0.99,
                                                       factor = 3))
})

test_that("the gap guard and the foundation's coverage give the same answer", {
  # One implementation, two callers. Sharing the 1.5x factor was not enough:
  # detection used to run on the whole series here and on the window's samples
  # there, so a gap straddling an edge was refused here and unseen there.
  tt <- c(0, .1, .2, .4, .45, .5, .6, .7, .8, .9, 1)   # a 0.2 s hole
  p <- PhysioExperiment::PhysioExperiment(
    assays = list(raw = matrix(seq_along(tt) * 1.0, length(tt), 1)),
    samplingRate = 10,
    rowData = S4Vectors::DataFrame(time_from_t0 = tt))
  m <- PhysioExperiment::MultiPhysioExperiment(streams = list(a = p),
                                               offsets = c(a = 0))

  windows <- list(c(0.25, 0.45),   # straddles the left edge
                  c(0.10, 0.30),   # straddles the right edge
                  c(0.25, 0.35),   # wholly inside the gap
                  c(0.40, 0.80),   # touches at 0.4 only
                  c(0.00, 0.20),   # touches at 0.2 only
                  c(0.50, 0.90))   # clear of it
  for (w in windows) {
    refused <- inherits(tryCatch(
      .reject_gap_in_window(tt, 10, "a", w[1], w[2]),
      error = function(e) e), "error")
    seen <- nrow(PhysioExperiment::streamCoverage(m, "a", w[1], w[2])$gaps) > 0L
    expect_identical(refused, seen,
                     info = sprintf("window [%g, %g]", w[1], w[2]))
  }
})
