# Tests for multi-channel coupling matrices (coupling-matrix.R)

# ---- coherenceMatrix ---------------------------------------------------------

test_that("coherenceMatrix returns correct dimensions", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- coherenceMatrix(mpe, modality_x = "EEG", modality_y = "EMG",
                            nperseg = 64L)

  expect_type(result, "list")
  expect_true("matrix" %in% names(result))
  expect_true("spectra" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true("channel_names_x" %in% names(result))
  expect_true("channel_names_y" %in% names(result))

  # 4 EEG channels x 2 EMG channels
  expect_equal(dim(result$matrix), c(4, 2))
  expect_equal(length(result$channel_names_x), 4)
  expect_equal(length(result$channel_names_y), 2)
})

test_that("coherenceMatrix values in [0, 1]", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- coherenceMatrix(mpe, modality_x = "EEG", modality_y = "EMG",
                            nperseg = 64L)

  expect_true(all(result$matrix >= 0 & result$matrix <= 1))
})

test_that("coherenceMatrix coupled channel has higher coherence", {
  pairs <- make_eeg_emg(n_sec = 4, sr_eeg = 200, sr_emg = 200,
                        cmc_strength = 0.9)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- coherenceMatrix(mpe, modality_x = "EEG", modality_y = "EMG",
                            nperseg = 128L)

  # Channel 1-1 is coupled; channel 2-2 is not
  expect_true(result$matrix[1, 1] > result$matrix[2, 2],
    info = sprintf("Coupled (1,1) = %.3f, uncoupled (2,2) = %.3f",
                   result$matrix[1, 1], result$matrix[2, 2]))
})

test_that("coherenceMatrix works with channel subset", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- coherenceMatrix(mpe, modality_x = "EEG", modality_y = "EMG",
                            channels_x = c(1L, 2L), channels_y = 1L,
                            nperseg = 64L)

  expect_equal(dim(result$matrix), c(2, 1))
})

test_that("coherenceMatrix with two PE objects", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)

  result <- coherenceMatrix(pairs$eeg, pairs$emg, nperseg = 64L)

  expect_type(result, "list")
  expect_equal(dim(result$matrix), c(4, 2))
})

test_that("coherenceMatrix preserves channel names from colData", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- coherenceMatrix(mpe, modality_x = "EEG", modality_y = "EMG",
                            nperseg = 64L)

  expect_true(all(grepl("^EEG", result$channel_names_x)))
  expect_true(all(grepl("^EMG", result$channel_names_y)))
})

# ---- couplingMatrix ----------------------------------------------------------

test_that("couplingMatrix returns correct dimensions with coherence", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- couplingMatrix(mpe, method = "coherence",
                           modality_x = "EEG", modality_y = "EMG",
                           nperseg = 64L)

  expect_type(result, "list")
  expect_true("matrix" %in% names(result))
  expect_true("method" %in% names(result))
  expect_equal(result$method, "coherence")
  expect_equal(dim(result$matrix), c(4, 2))
})

test_that("couplingMatrix works with crosscorrelation method", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- couplingMatrix(mpe, method = "crosscorrelation",
                           modality_x = "EEG", modality_y = "EMG")

  expect_type(result, "list")
  expect_equal(dim(result$matrix), c(4, 2))
  expect_true(all(result$matrix >= 0))
})

test_that("couplingMatrix compatible with plotCouplingMatrix", {
  skip_if_not_installed("ggplot2")

  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- couplingMatrix(mpe, method = "coherence",
                           modality_x = "EEG", modality_y = "EMG",
                           nperseg = 64L)

  p <- plotCouplingMatrix(result$matrix)
  expect_s3_class(p, "ggplot")
})

test_that("couplingMatrix with channel subset", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- couplingMatrix(mpe, method = "coherence",
                           modality_x = "EEG", modality_y = "EMG",
                           channels_x = 1L, channels_y = c(1L, 2L),
                           nperseg = 64L)

  expect_equal(dim(result$matrix), c(1, 2))
})

# ---- .resolve_pe_pair --------------------------------------------------------

test_that(".resolve_pe_pair works with MPE", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  result <- PhysioCrossModal:::.resolve_pe_pair(mpe, NULL, "EEG", "EMG")
  expect_s4_class(result$pe_x, "PhysioExperiment")
  expect_s4_class(result$pe_y, "PhysioExperiment")
})

test_that(".resolve_pe_pair works with two PEs", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)

  result <- PhysioCrossModal:::.resolve_pe_pair(pairs$eeg, pairs$emg,
                                                 NULL, NULL)
  expect_s4_class(result$pe_x, "PhysioExperiment")
  expect_s4_class(result$pe_y, "PhysioExperiment")
})

test_that(".resolve_pe_pair errors on invalid input", {
  expect_error(
    PhysioCrossModal:::.resolve_pe_pair("not_pe", NULL, NULL, NULL),
    "must be"
  )
})

test_that(".resolve_pe_pair requires modalities for MPE", {
  pairs <- make_eeg_emg(n_sec = 2, sr_eeg = 200, sr_emg = 200)
  mpe <- MultiPhysioExperiment(EEG = pairs$eeg, EMG = pairs$emg)

  expect_error(
    PhysioCrossModal:::.resolve_pe_pair(mpe, NULL, NULL, NULL),
    "modality_x and modality_y"
  )
})
