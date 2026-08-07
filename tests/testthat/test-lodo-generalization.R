# Tests for LODO generalization benchmark

test_that("lodoGeneralization works for gaussian outcomes", {
  set.seed(1001)
  n <- 180
  df <- data.frame(
    site_id = rep(c("A", "B", "C"), each = 60),
    f1 = rnorm(n),
    f2 = rnorm(n),
    stringsAsFactors = FALSE
  )
  df$y <- 1.2 * df$f1 - 0.7 * df$f2 + rnorm(n, sd = 0.3)

  out <- lodoGeneralization(
    data = df,
    outcome = "y",
    site = "site_id",
    family = "gaussian"
  )

  expect_type(out, "list")
  expect_true(all(c("fold_metrics", "predictions", "aggregate", "settings") %in%
                  names(out)))
  expect_equal(nrow(out$fold_metrics), 3)
  expect_true(all(out$fold_metrics$status == "ok"))
  expect_true(all(c("rmse", "mae", "r2") %in% names(out$fold_metrics)))
  expect_true(is.finite(out$aggregate$rmse))
})

test_that("lodoGeneralization works for binomial outcomes", {
  set.seed(1002)
  n <- 210
  df <- data.frame(
    site_id = rep(c("A", "B", "C"), each = 70),
    f1 = rnorm(n),
    f2 = rnorm(n),
    stringsAsFactors = FALSE
  )

  lin <- 1.1 * df$f1 - 0.8 * df$f2 + ifelse(df$site_id == "C", 0.2, 0)
  prob <- 1 / (1 + exp(-lin))
  df$improved <- ifelse(stats::runif(n) < prob, "yes", "no")

  out <- lodoGeneralization(
    data = df,
    outcome = "improved",
    site = "site_id",
    family = "binomial",
    positive_class = "yes"
  )

  expect_equal(nrow(out$fold_metrics), 3)
  expect_true(all(out$fold_metrics$status == "ok"))
  expect_true(all(c("accuracy", "balanced_accuracy", "auc") %in%
                  names(out$fold_metrics)))
  expect_true(all(out$fold_metrics$auc >= 0 & out$fold_metrics$auc <= 1))
})

test_that("lodoGeneralization validates minimum requirements", {
  df <- data.frame(
    site_id = rep("A", 20),
    f1 = rnorm(20),
    y = rnorm(20),
    stringsAsFactors = FALSE
  )

  expect_error(
    lodoGeneralization(df, outcome = "y", site = "site_id"),
    "at least 2 sites"
  )

  expect_error(
    lodoGeneralization(df, outcome = "y", site = "site_id", features = "missing"),
    "not present"
  )
})
