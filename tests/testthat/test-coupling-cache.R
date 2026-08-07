make_mpe <- function() {
  eeg <- PhysioExperiment(assays = list(raw = matrix(rnorm(1000), 500, 2)),
                          samplingRate = 250)
  emg <- PhysioExperiment(assays = list(raw = matrix(rnorm(1000), 500, 2)),
                          samplingRate = 250)
  MultiPhysioExperiment(EEG = eeg, EMG = emg)
}

test_that("couplingResults getter/setter round-trips", {
  mpe <- make_mpe()
  expect_length(couplingResults(mpe), 0L)
  couplingResults(mpe) <- list(a = 1)
  expect_identical(couplingResults(mpe)$a, 1)
})

test_that("couplingAnalysisCached memoises by argument digest (cache hit, no recompute)", {
  set.seed(1)
  mpe <- make_mpe()

  r1 <- couplingAnalysisCached(mpe, "coherence", "EEG", "EMG")
  expect_false(r1$cached)
  expect_equal(length(couplingResults(r1$mpe)), 1L)   # cache grew once

  r2 <- couplingAnalysisCached(r1$mpe, "coherence", "EEG", "EMG")
  expect_true(r2$cached)                              # served from cache
  expect_identical(r1$result, r2$result)             # byte-identical
  expect_equal(length(couplingResults(r2$mpe)), 1L)   # no new entry

  # different arguments -> distinct key -> new entry
  r3 <- couplingAnalysisCached(r2$mpe, "coherence", "EEG", "EMG", channels_x = 2L)
  expect_false(r3$cached)
  expect_equal(length(couplingResults(r3$mpe)), 2L)
})

test_that("cache = FALSE forces recompute but still stores the result", {
  set.seed(2)
  mpe <- make_mpe()
  r1 <- couplingAnalysisCached(mpe, "coherence", "EEG", "EMG")
  r2 <- couplingAnalysisCached(r1$mpe, "coherence", "EEG", "EMG", cache = FALSE)
  expect_false(r2$cached)                              # not served from cache
  expect_equal(length(couplingResults(r2$mpe)), 1L)    # key overwritten, not duplicated
})
