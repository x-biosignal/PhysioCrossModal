# Re-source data helpers after load_all() so they shadow same-named examples.
source(testthat::test_path("helper-data.R"), local = TRUE)
