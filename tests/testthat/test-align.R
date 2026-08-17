library(testthat)
library(PhysioCrossModal)

# ---- helpers ----------------------------------------------------------------

make_pe <- function(n_time = 1000, n_ch = 4, sr = 1000, seed = 42) {
  set.seed(seed)
  data <- matrix(rnorm(n_time * n_ch), nrow = n_time, ncol = n_ch)
  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = paste0("Ch", seq_len(n_ch)),
      type  = rep("EEG", n_ch)
    ),
    samplingRate = sr
  )
}

# ---- alignToRate tests ------------------------------------------------------

test_that("alignToRate resamples 1000Hz PE to 500Hz correctly", {
  pe <- make_pe(n_time = 1000, n_ch = 4, sr = 1000)
  pe_out <- alignToRate(pe, target_rate = 500)

  expect_s4_class(pe_out, "PhysioExperiment")
  expect_equal(samplingRate(pe_out), 500)
  # 1 s of data at 1000Hz -> 1 s at 500Hz = 500 samples
  expect_equal(nrow(SummarizedExperiment::assay(pe_out)), 500)
  # Channels preserved

  expect_equal(ncol(SummarizedExperiment::assay(pe_out)), 4)
})

test_that("alignToRate returns input unchanged when rate already matches", {
  pe <- make_pe(n_time = 500, n_ch = 2, sr = 500)
  pe_out <- alignToRate(pe, target_rate = 500)

  expect_identical(pe_out, pe)
})

test_that("alignToRate works with 3D assay", {
  set.seed(123)
  arr <- array(rnorm(1000 * 4 * 3), dim = c(1000, 4, 3))
  pe <- PhysioExperiment(
    assays = list(raw = arr),
    samplingRate = 1000
  )
  pe_out <- alignToRate(pe, target_rate = 500)

  expect_equal(samplingRate(pe_out), 500)
  out_dims <- dim(SummarizedExperiment::assay(pe_out))
  expect_equal(out_dims[1], 500)
  expect_equal(out_dims[2], 4)
  expect_equal(out_dims[3], 3)
})

test_that("alignToRate validates inputs", {
  pe <- make_pe()
  expect_error(alignToRate(pe, target_rate = -1), "positive")
  expect_error(alignToRate(pe, target_rate = "abc"), "positive")
  expect_error(alignToRate("not_pe", target_rate = 500), "PhysioExperiment")
})

# ---- alignSignals tests ----------------------------------------------------

test_that("alignSignals with lowest_rate: 500Hz + 1000Hz -> both 500Hz", {
  pe500  <- make_pe(n_time = 500,  n_ch = 2, sr = 500)
  pe1000 <- make_pe(n_time = 1000, n_ch = 3, sr = 1000)

  mpe <- alignSignals(EEG = pe500, EMG = pe1000, method = "lowest_rate")

  expect_s4_class(mpe, "MultiPhysioExperiment")
  rates <- samplingRates(mpe)
  expect_equal(unname(rates), c(500, 500))
})

test_that("alignSignals with target_rate: 500Hz + 1000Hz -> both 250Hz", {
  pe500  <- make_pe(n_time = 500,  n_ch = 2, sr = 500)
  pe1000 <- make_pe(n_time = 1000, n_ch = 3, sr = 1000)

  mpe <- alignSignals(EEG = pe500, EMG = pe1000,
                       method = "resample", target_rate = 250)

  expect_s4_class(mpe, "MultiPhysioExperiment")
  rates <- samplingRates(mpe)
  expect_equal(unname(rates), c(250, 250))

  # Both should have duration * 250 samples
  eeg_assay <- SummarizedExperiment::assay(experiments(mpe)[["EEG"]])
  emg_assay <- SummarizedExperiment::assay(experiments(mpe)[["EMG"]])
  expect_equal(nrow(eeg_assay), 250)
  expect_equal(nrow(emg_assay), 250)
})

test_that("alignSignals same rate returns MPE without resampling", {
  pe1 <- make_pe(n_time = 500, n_ch = 2, sr = 500, seed = 1)
  pe2 <- make_pe(n_time = 500, n_ch = 3, sr = 500, seed = 2)

  mpe <- alignSignals(A = pe1, B = pe2)

  expect_s4_class(mpe, "MultiPhysioExperiment")
  rates <- samplingRates(mpe)
  expect_equal(unname(rates), c(500, 500))

  # Data should be untouched
  expect_equal(
    SummarizedExperiment::assay(experiments(mpe)[["A"]]),
    SummarizedExperiment::assay(pe1)
  )
  expect_equal(
    SummarizedExperiment::assay(experiments(mpe)[["B"]]),
    SummarizedExperiment::assay(pe2)
  )
})

test_that("alignSignals with common_rate picks highest rate", {
  pe500  <- make_pe(n_time = 500,  n_ch = 2, sr = 500)
  pe1000 <- make_pe(n_time = 1000, n_ch = 3, sr = 1000)

  mpe <- alignSignals(EEG = pe500, EMG = pe1000, method = "common_rate")

  rates <- samplingRates(mpe)
  expect_equal(unname(rates), c(1000, 1000))
})

test_that("alignSignals requires target_rate for method='resample'", {
  pe <- make_pe(n_time = 500, n_ch = 2, sr = 500)
  expect_error(
    alignSignals(A = pe, method = "resample"),
    "target_rate is required"
  )
})

test_that("alignSignals requires named arguments", {
  pe <- make_pe()
  expect_error(alignSignals(pe), "named")
})

# ---- mergePhysio tests ------------------------------------------------------

test_that("mergePhysio combines two 4-ch PEs into 8-ch PE with prefixed labels", {
  pe1 <- make_pe(n_time = 100, n_ch = 4, sr = 500, seed = 1)
  pe2 <- make_pe(n_time = 100, n_ch = 4, sr = 500, seed = 2)

  merged <- mergePhysio(pe1, pe2)

  expect_s4_class(merged, "PhysioExperiment")
  expect_equal(ncol(SummarizedExperiment::assay(merged)), 8)
  expect_equal(nrow(SummarizedExperiment::assay(merged)), 100)
  expect_equal(samplingRate(merged), 500)

  labels <- SummarizedExperiment::colData(merged)$label
  expect_equal(length(labels), 8)
  # First four should have x_ prefix, last four y_ prefix
  expect_true(all(grepl("^x_", labels[1:4])))
  expect_true(all(grepl("^y_", labels[5:8])))
})

test_that("mergePhysio errors when sampling rates differ", {
  pe1 <- make_pe(n_time = 100, n_ch = 2, sr = 500)
  pe2 <- make_pe(n_time = 100, n_ch = 2, sr = 1000)

  expect_error(mergePhysio(pe1, pe2), "Sampling rates must match")
})

test_that("mergePhysio errors when time points differ", {
  pe1 <- make_pe(n_time = 100, n_ch = 2, sr = 500)
  pe2 <- make_pe(n_time = 200, n_ch = 2, sr = 500)

  expect_error(mergePhysio(pe1, pe2), "time points must match")
})

# ---- .resample_fft tests -----------------------------------------------------

test_that(".resample_fft odd downsample produces correct length without warning", {
  x <- sin(2 * pi * 5 * seq(0, 1, length.out = 100))
  # Downsample from 100 to 77 (odd target)
  result <- PhysioCrossModal:::.resample_fft(x, 77)
  expect_length(result, 77)
})

test_that(".resample_fft even downsample produces correct length without warning", {
  x <- sin(2 * pi * 5 * seq(0, 1, length.out = 100))
  # Downsample from 100 to 60 (even target)
  result <- PhysioCrossModal:::.resample_fft(x, 60)
  expect_length(result, 60)
})

test_that(".resample_fft roundtrip upsample-then-downsample recovers signal", {
  set.seed(42)
  n_orig <- 100
  x <- sin(2 * pi * 3 * seq(0, 1, length.out = n_orig))

  # Upsample to 200, then back to 100
  upsampled <- PhysioCrossModal:::.resample_fft(x, 200)
  roundtrip <- PhysioCrossModal:::.resample_fft(upsampled, n_orig)

  expect_length(roundtrip, n_orig)
  # Should recover original signal closely
  expect_equal(roundtrip, x, tolerance = 0.05)
})

test_that(".resample_fft identity returns input unchanged", {
  x <- rnorm(50)
  result <- PhysioCrossModal:::.resample_fft(x, 50)
  expect_equal(result, x)
})

test_that("mergePhysio respects custom prefix", {
  pe1 <- make_pe(n_time = 100, n_ch = 2, sr = 500, seed = 1)
  pe2 <- make_pe(n_time = 100, n_ch = 2, sr = 500, seed = 2)

  merged <- mergePhysio(pe1, pe2, prefix = c("eeg_", "emg_"))

  labels <- SummarizedExperiment::colData(merged)$label
  expect_true(all(grepl("^eeg_", labels[1:2])))
  expect_true(all(grepl("^emg_", labels[3:4])))
})
