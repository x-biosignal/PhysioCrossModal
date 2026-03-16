library(testthat)
library(PhysioCrossModal)

# ---- Helper: create PhysioExperiment objects for tests --------------------

make_test_pe <- function(n_time = 500, n_ch = 4, sr = 250,
                          ch_prefix = "Ch") {
  data <- matrix(rnorm(n_time * n_ch), nrow = n_time, ncol = n_ch)
  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = paste0(ch_prefix, seq_len(n_ch)),
      type = rep(ch_prefix, n_ch)
    ),
    samplingRate = sr
  )
}

# ---- Constructor tests ---------------------------------------------------

test_that("Constructor works with named PE objects at different sampling rates", {
  pe_eeg <- make_test_pe(n_time = 500, n_ch = 4, sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(n_time = 1000, n_ch = 2, sr = 1000, ch_prefix = "EMG")

  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  expect_s4_class(mpe, "MultiPhysioExperiment")
  expect_equal(nModalities(mpe), 2L)
  expect_equal(modalities(mpe), c("EEG", "EMG"))
})

test_that("Constructor works with experiments list argument", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_ecg <- make_test_pe(n_ch = 1, sr = 500, ch_prefix = "ECG")

  mpe <- MultiPhysioExperiment(experiments = list(EEG = pe_eeg, ECG = pe_ecg))

  expect_s4_class(mpe, "MultiPhysioExperiment")
  expect_equal(nModalities(mpe), 2L)
})

test_that("Constructor works with single PE object", {
  pe <- make_test_pe(sr = 250, ch_prefix = "EEG")
  mpe <- MultiPhysioExperiment(EEG = pe)

  expect_s4_class(mpe, "MultiPhysioExperiment")
  expect_equal(nModalities(mpe), 1L)
})

test_that("Constructor accepts optional alignment DataFrame", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")

  aln <- S4Vectors::DataFrame(
    modality = c("EEG", "EMG"),
    samplingRate = c(250, 1000),
    offset = c(0, 0.5)
  )

  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg, alignment = aln)

  expect_s4_class(mpe, "MultiPhysioExperiment")
  expect_equal(nrow(alignment(mpe)), 2L)
  expect_equal(alignment(mpe)$offset, c(0, 0.5))
})

# ---- Validation tests ----------------------------------------------------

test_that("Constructor rejects unnamed experiments", {
  pe1 <- make_test_pe(sr = 250)
  pe2 <- make_test_pe(sr = 500)

  expect_error(
    MultiPhysioExperiment(pe1, pe2),
    "all experiments must be named"
  )
})

test_that("Constructor rejects non-PhysioExperiment objects", {
  expect_error(
    MultiPhysioExperiment(EEG = data.frame(x = 1:10)),
    "all experiments must be PhysioExperiment objects"
  )
})

test_that("Constructor rejects empty experiments list", {
  expect_error(
    MultiPhysioExperiment(),
    "experiments must contain at least one PhysioExperiment"
  )
})

test_that("Constructor rejects partially unnamed experiments", {
  pe1 <- make_test_pe(sr = 250)
  pe2 <- make_test_pe(sr = 500)

  expect_error(
    MultiPhysioExperiment(EEG = pe1, pe2),
    "all experiments must be named"
  )
})

# ---- Accessor tests ------------------------------------------------------

test_that("experiments() returns named list of PhysioExperiment objects", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  exps <- experiments(mpe)

  expect_type(exps, "list")
  expect_named(exps, c("EEG", "EMG"))
  expect_s4_class(exps[["EEG"]], "PhysioExperiment")
  expect_s4_class(exps[["EMG"]], "PhysioExperiment")
})

test_that("experiments<-() replacement works and validates", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg)

  experiments(mpe) <- list(EEG = pe_eeg, EMG = pe_emg)

  expect_equal(nModalities(mpe), 2L)

  # Replacement with invalid data should error
  expect_error(
    experiments(mpe) <- list(EEG = data.frame(x = 1)),
    "all experiments must be PhysioExperiment objects"
  )
})

test_that("modalities() returns correct character vector", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  pe_ecg <- make_test_pe(n_ch = 1, sr = 500, ch_prefix = "ECG")

  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg, ECG = pe_ecg)

  mods <- modalities(mpe)

  expect_type(mods, "character")
  expect_equal(length(mods), 3L)
  expect_true("EEG" %in% mods)
  expect_true("EMG" %in% mods)
  expect_true("ECG" %in% mods)
})

test_that("samplingRates() returns correct named numeric vector", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")

  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  srs <- samplingRates(mpe)

  expect_type(srs, "double")
  expect_named(srs, c("EEG", "EMG"))
  expect_equal(srs[["EEG"]], 250)
  expect_equal(srs[["EMG"]], 1000)
})

test_that("nModalities() returns correct integer", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg)
  expect_equal(nModalities(mpe), 1L)

  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  mpe2 <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)
  expect_equal(nModalities(mpe2), 2L)
})

test_that("alignment() and alignment<-() work correctly", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  aln <- alignment(mpe)
  expect_s4_class(aln, "DataFrame")

  # Set new alignment
  new_aln <- S4Vectors::DataFrame(
    modality = c("EEG", "EMG"),
    samplingRate = c(250, 1000),
    offset = c(0, 1.0)
  )
  alignment(mpe) <- new_aln
  expect_equal(alignment(mpe)$offset, c(0, 1.0))
})

# ---- show() test ---------------------------------------------------------

test_that("show() outputs 'MultiPhysioExperiment'", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  output <- capture.output(show(mpe))

  expect_true(any(grepl("MultiPhysioExperiment", output)))
  expect_true(any(grepl("EEG", output)))
  expect_true(any(grepl("EMG", output)))
  expect_true(any(grepl("250", output)))
  expect_true(any(grepl("1000", output)))
})

# ---- [[ subsetting -------------------------------------------------------

test_that("[[ extracts a single PhysioExperiment by name", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  extracted <- mpe[["EEG"]]
  expect_s4_class(extracted, "PhysioExperiment")
  expect_equal(samplingRate(extracted), 250)
})

test_that("[[ extracts a single PhysioExperiment by index", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  extracted <- mpe[[1]]
  expect_s4_class(extracted, "PhysioExperiment")
  expect_equal(samplingRate(extracted), 250)
})

# ---- Default alignment ---------------------------------------------------

test_that("default alignment is built from experiments", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  aln <- alignment(mpe)

  expect_s4_class(aln, "DataFrame")
  expect_equal(nrow(aln), 2L)
  expect_true("modality" %in% names(aln))
  expect_true("samplingRate" %in% names(aln))
  expect_equal(aln$modality, c("EEG", "EMG"))
  expect_equal(unname(aln$samplingRate), c(250, 1000))
})

# ---- length() method -----------------------------------------------------

test_that("length() returns number of modalities", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  mpe1 <- MultiPhysioExperiment(EEG = pe_eeg)
  expect_equal(length(mpe1), 1L)

  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  pe_ecg <- make_test_pe(n_ch = 1, sr = 500, ch_prefix = "ECG")
  mpe3 <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg, ECG = pe_ecg)
  expect_equal(length(mpe3), 3L)
  expect_equal(length(mpe3), nModalities(mpe3))
})

# ---- names() method ------------------------------------------------------

test_that("names() returns modality names", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  nms <- names(mpe)
  expect_type(nms, "character")
  expect_equal(nms, c("EEG", "EMG"))
  expect_equal(nms, modalities(mpe))
})

# ---- [ modality subsetting -----------------------------------------------

test_that("mpe[, 'EEG'] subsets by single modality name", {
  pe_eeg <- make_test_pe(n_time = 500, n_ch = 4, sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(n_time = 1000, n_ch = 2, sr = 1000, ch_prefix = "EMG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  mpe_eeg <- mpe[, "EEG"]

  expect_s4_class(mpe_eeg, "MultiPhysioExperiment")
  expect_equal(nModalities(mpe_eeg), 1L)
  expect_equal(modalities(mpe_eeg), "EEG")
  expect_equal(samplingRates(mpe_eeg)[["EEG"]], 250)
})

test_that("mpe[, c('EEG', 'EMG')] subsets by multiple modality names", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(sr = 1000, ch_prefix = "EMG")
  pe_ecg <- make_test_pe(n_ch = 1, sr = 500, ch_prefix = "ECG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg, ECG = pe_ecg)

  mpe_sub <- mpe[, c("EEG", "EMG")]

  expect_s4_class(mpe_sub, "MultiPhysioExperiment")
  expect_equal(nModalities(mpe_sub), 2L)
  expect_equal(modalities(mpe_sub), c("EEG", "EMG"))
})

test_that("mpe[, 'invalid'] errors for unknown modality names", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg)

  expect_error(mpe[, "NONEXISTENT"], "Unknown modality name")
  expect_error(mpe[, c("EEG", "FAKE")], "Unknown modality name")
})

# ---- [ time-range subsetting ---------------------------------------------

test_that("mpe[c(1.0, 5.0), ] subsets by time range respecting sampling rates", {
  # EEG: 500 Hz, 5 seconds = 2500 samples
  pe_eeg <- make_test_pe(n_time = 2500, n_ch = 4, sr = 500, ch_prefix = "EEG")
  # EMG: 1000 Hz, 2.5 seconds = 2500 samples
  pe_emg <- make_test_pe(n_time = 2500, n_ch = 2, sr = 1000, ch_prefix = "EMG")

  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  # Subset time range [1.0, 3.0) seconds
  mpe_sub <- mpe[c(1.0, 3.0), ]

  expect_s4_class(mpe_sub, "MultiPhysioExperiment")
  expect_equal(nModalities(mpe_sub), 2L)

  # EEG@500Hz: start = floor(1.0*500)+1 = 501, end = floor(3.0*500)+1 = 1501
  # -> 1501 - 501 + 1 = 1001 samples
  eeg_sub <- experiments(mpe_sub)[["EEG"]]
  expect_equal(length(eeg_sub), 1001L)

  # EMG@1000Hz: start = floor(1.0*1000)+1 = 1001, end = floor(3.0*1000)+1 = 2500 (clamped)
  # -> min(2500, 3001) = 2500, so end = 2500
  # Wait, floor(3.0*1000)+1 = 3001, but n_samples = 2500, so clamped to 2500
  # -> 2500 - 1001 + 1 = 1500 samples
  emg_sub <- experiments(mpe_sub)[["EMG"]]
  expect_equal(length(emg_sub), 1500L)
})

test_that("time-range subsetting with exact sample-count verification", {
  # Create data with known dimensions for easy math
  # EEG: 250 Hz, 2000 samples = 8 seconds
  pe_eeg <- make_test_pe(n_time = 2000, n_ch = 2, sr = 250, ch_prefix = "EEG")
  # EMG: 1000 Hz, 8000 samples = 8 seconds
  pe_emg <- make_test_pe(n_time = 8000, n_ch = 2, sr = 1000, ch_prefix = "EMG")

  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  # Subset [2.0, 4.0]
  mpe_sub <- mpe[c(2.0, 4.0), ]

  # EEG@250Hz: start = floor(2.0*250)+1 = 501, end = floor(4.0*250)+1 = 1001
  # -> 1001 - 501 + 1 = 501 samples
  expect_equal(length(experiments(mpe_sub)[["EEG"]]), 501L)

  # EMG@1000Hz: start = floor(2.0*1000)+1 = 2001, end = floor(4.0*1000)+1 = 4001
  # -> 4001 - 2001 + 1 = 2001 samples
  expect_equal(length(experiments(mpe_sub)[["EMG"]]), 2001L)
})

# ---- [ combined subsetting -----------------------------------------------

test_that("mpe[c(1.0, 3.0), 'EEG'] subsets by both time range and modality", {
  pe_eeg <- make_test_pe(n_time = 2000, n_ch = 4, sr = 250, ch_prefix = "EEG")
  pe_emg <- make_test_pe(n_time = 8000, n_ch = 2, sr = 1000, ch_prefix = "EMG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

  mpe_sub <- mpe[c(1.0, 3.0), "EEG"]

  expect_s4_class(mpe_sub, "MultiPhysioExperiment")
  expect_equal(nModalities(mpe_sub), 1L)
  expect_equal(modalities(mpe_sub), "EEG")

  # EEG@250Hz: start = floor(1.0*250)+1 = 251, end = floor(3.0*250)+1 = 751
  # -> 751 - 251 + 1 = 501 samples
  expect_equal(length(experiments(mpe_sub)[["EEG"]]), 501L)
})

# ---- [ input validation --------------------------------------------------

test_that("time-range subsetting errors on invalid i", {
  pe_eeg <- make_test_pe(sr = 250, ch_prefix = "EEG")
  mpe <- MultiPhysioExperiment(EEG = pe_eeg)

  # i must be length 2

  expect_error(mpe[1.0, ], "numeric vector of length 2")

  # i must be numeric
  expect_error(mpe[c("a", "b"), ], "numeric vector of length 2")

  # tmin must be less than tmax
  expect_error(mpe[c(5.0, 1.0), ], "tmin must be less than tmax")
})
