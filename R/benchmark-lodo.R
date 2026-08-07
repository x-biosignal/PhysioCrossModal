# External validity benchmarking utilities
#
# Implements leave-one-site-out (LODO) benchmarking to quantify
# cross-site generalization in rehabilitation sensing studies.

#' Leave-one-site-out generalization benchmark
#'
#' Evaluates model transportability by training on all-but-one sites and
#' testing on the held-out site (LODO). Supports continuous outcomes
#' (`family = "gaussian"`) and binary outcomes (`family = "binomial"`).
#'
#' @param data Data frame containing outcome, site label, and features.
#' @param outcome Outcome column name.
#' @param site Site/facility column name used for LODO splitting.
#' @param features Feature column names. If `NULL`, numeric columns excluding
#'   `outcome` and `site` are used.
#' @param family Modeling family: `"gaussian"` or `"binomial"`.
#' @param positive_class Positive class label for binomial metrics.
#'   If `NULL`, the second factor level is used.
#' @param threshold Classification threshold for binomial predictions.
#' @param min_train_rows Minimum training rows required per fold.
#' @param scale_features Logical; z-score features using training statistics.
#'
#' @return A list with:
#' \describe{
#'   \item{fold_metrics}{Per-site metrics.}
#'   \item{predictions}{Row-level held-out predictions.}
#'   \item{aggregate}{Mean metrics across folds.}
#'   \item{settings}{Benchmark settings and feature list.}
#' }
#' @seealso [couplingAnalysis()], [surrogateTest()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 120
#' df <- data.frame(
#'   site_id = rep(c("A", "B", "C"), each = 40),
#'   f1 = rnorm(n),
#'   f2 = rnorm(n)
#' )
#' df$y <- 0.8 * df$f1 - 0.4 * df$f2 + rnorm(n, sd = 0.2)
#'
#' res <- lodoGeneralization(
#'   data = df,
#'   outcome = "y",
#'   site = "site_id",
#'   family = "gaussian"
#' )
#' head(res$fold_metrics)
lodoGeneralization <- function(data,
                               outcome,
                               site,
                               features = NULL,
                               family = c("gaussian", "binomial"),
                               positive_class = NULL,
                               threshold = 0.5,
                               min_train_rows = 20L,
                               scale_features = TRUE) {
  family <- match.arg(family)

  if (!is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }
  if (!all(c(outcome, site) %in% names(data))) {
    stop("outcome/site columns not found in data", call. = FALSE)
  }

  if (is.null(features)) {
    numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
    features <- setdiff(numeric_cols, c(outcome, site))
  }
  if (length(features) == 0L) {
    stop("No features available for modeling", call. = FALSE)
  }
  if (!all(features %in% names(data))) {
    stop("Some features are not present in data", call. = FALSE)
  }

  keep_cols <- unique(c(site, outcome, features))
  df <- data[, keep_cols, drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]

  sites <- unique(as.character(df[[site]]))
  if (length(sites) < 2L) {
    stop("LODO requires at least 2 sites", call. = FALSE)
  }

  min_train_rows <- as.integer(min_train_rows)
  fold_metrics <- vector("list", length(sites))
  pred_rows <- vector("list", length(sites))

  for (i in seq_along(sites)) {
    holdout <- sites[[i]]
    idx_test <- df[[site]] == holdout
    idx_train <- !idx_test

    train <- df[idx_train, , drop = FALSE]
    test <- df[idx_test, , drop = FALSE]

    if (nrow(train) < min_train_rows || nrow(test) == 0L) {
      if (family == "gaussian") {
        fold_metrics[[i]] <- data.frame(
          holdout_site = holdout,
          n_train = nrow(train),
          n_test = nrow(test),
          status = "skipped",
          rmse = NA_real_,
          mae = NA_real_,
          r2 = NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        fold_metrics[[i]] <- data.frame(
          holdout_site = holdout,
          n_train = nrow(train),
          n_test = nrow(test),
          status = "skipped",
          accuracy = NA_real_,
          sensitivity = NA_real_,
          specificity = NA_real_,
          balanced_accuracy = NA_real_,
          auc = NA_real_,
          stringsAsFactors = FALSE
        )
      }
      pred_rows[[i]] <- data.frame()
      next
    }

    x_train <- as.matrix(train[, features, drop = FALSE])
    x_test <- as.matrix(test[, features, drop = FALSE])

    if (isTRUE(scale_features)) {
      mu <- colMeans(x_train)
      sigma <- apply(x_train, 2, stats::sd)
      sigma[sigma == 0] <- 1
      x_train <- sweep(sweep(x_train, 2, mu, FUN = "-"), 2, sigma, FUN = "/")
      x_test <- sweep(sweep(x_test, 2, mu, FUN = "-"), 2, sigma, FUN = "/")
    }

    dat_train <- as.data.frame(x_train, stringsAsFactors = FALSE)
    dat_test <- as.data.frame(x_test, stringsAsFactors = FALSE)

    form <- stats::as.formula(paste("y ~", paste(features, collapse = " + ")))

    if (family == "gaussian") {
      dat_train$y <- as.numeric(train[[outcome]])
      dat_test$y <- as.numeric(test[[outcome]])

      fit <- stats::lm(form, data = dat_train)
      pred <- as.numeric(stats::predict(fit, newdata = dat_test))
      obs <- dat_test$y

      rmse <- sqrt(mean((pred - obs)^2))
      mae <- mean(abs(pred - obs))
      r2 <- if (stats::sd(obs) > 0) stats::cor(pred, obs)^2 else NA_real_

      fold_metrics[[i]] <- data.frame(
        holdout_site = holdout,
        n_train = nrow(train),
        n_test = nrow(test),
        status = "ok",
        rmse = rmse,
        mae = mae,
        r2 = r2,
        stringsAsFactors = FALSE
      )

      pred_rows[[i]] <- data.frame(
        holdout_site = holdout,
        observed = obs,
        predicted = pred,
        stringsAsFactors = FALSE
      )
    } else {
      y_train <- train[[outcome]]
      y_test <- test[[outcome]]
      if (is.logical(y_train)) y_train <- as.integer(y_train)
      if (is.logical(y_test)) y_test <- as.integer(y_test)

      y_train <- as.factor(y_train)
      y_test <- as.factor(y_test)

      if (nlevels(y_train) != 2L) {
        stop("Binomial outcome must have exactly 2 classes", call. = FALSE)
      }

      if (is.null(positive_class)) {
        positive <- levels(y_train)[2]
      } else {
        positive <- as.character(positive_class)
        if (!positive %in% levels(y_train)) {
          stop("positive_class not found in training data levels", call. = FALSE)
        }
      }

      dat_train$y <- stats::relevel(y_train, ref = setdiff(levels(y_train), positive)[1])
      dat_test$y <- factor(y_test, levels = levels(dat_train$y))

      fit <- stats::glm(form, data = dat_train, family = stats::binomial())
      prob <- as.numeric(stats::predict(fit, newdata = dat_test, type = "response"))

      pred_class <- ifelse(prob >= threshold, positive, setdiff(levels(dat_train$y), positive)[1])
      obs_chr <- as.character(dat_test$y)

      tp <- sum(pred_class == positive & obs_chr == positive, na.rm = TRUE)
      tn <- sum(pred_class != positive & obs_chr != positive, na.rm = TRUE)
      fp <- sum(pred_class == positive & obs_chr != positive, na.rm = TRUE)
      fn <- sum(pred_class != positive & obs_chr == positive, na.rm = TRUE)

      accuracy <- (tp + tn) / max(1, tp + tn + fp + fn)
      sensitivity <- tp / max(1, tp + fn)
      specificity <- tn / max(1, tn + fp)
      bal_acc <- 0.5 * (sensitivity + specificity)
      auc <- .binary_auc(obs_chr == positive, prob)

      fold_metrics[[i]] <- data.frame(
        holdout_site = holdout,
        n_train = nrow(train),
        n_test = nrow(test),
        status = "ok",
        accuracy = accuracy,
        sensitivity = sensitivity,
        specificity = specificity,
        balanced_accuracy = bal_acc,
        auc = auc,
        stringsAsFactors = FALSE
      )

      pred_rows[[i]] <- data.frame(
        holdout_site = holdout,
        observed = obs_chr,
        predicted = pred_class,
        predicted_prob = prob,
        positive_class = positive,
        stringsAsFactors = FALSE
      )
    }
  }

  fold_metrics_df <- do.call(rbind, fold_metrics)
  pred_df <- do.call(rbind, pred_rows)

  ok_rows <- fold_metrics_df$status == "ok"
  metric_cols <- setdiff(
    names(fold_metrics_df),
    c("holdout_site", "n_train", "n_test", "status")
  )
  agg <- list()
  for (col in metric_cols) {
    agg[[col]] <- if (any(ok_rows)) {
      mean(fold_metrics_df[[col]][ok_rows], na.rm = TRUE)
    } else {
      NA_real_
    }
  }
  agg$n_folds_ok <- sum(ok_rows)
  agg$n_folds_total <- nrow(fold_metrics_df)

  list(
    fold_metrics = fold_metrics_df,
    predictions = pred_df,
    aggregate = agg,
    settings = list(
      outcome = outcome,
      site = site,
      features = features,
      family = family,
      threshold = threshold,
      min_train_rows = min_train_rows,
      scale_features = scale_features
    )
  )
}


#' AUC for binary classification
#'
#' Rank-based AUC estimate equivalent to the Mann-Whitney statistic.
#'
#' @param truth Logical vector where TRUE indicates positive class.
#' @param score Numeric prediction scores.
#' @return Numeric AUC in [0,1], or NA when undefined.
#' @noRd
.binary_auc <- function(truth, score) {
  truth <- as.logical(truth)
  if (length(truth) != length(score)) return(NA_real_)
  keep <- !is.na(truth) & !is.na(score)
  truth <- truth[keep]
  score <- score[keep]
  if (!length(truth)) return(NA_real_)
  if (all(!truth) || all(truth)) return(NA_real_)

  ranks <- rank(score, ties.method = "average")
  n_pos <- sum(truth)
  n_neg <- sum(!truth)
  (sum(ranks[truth]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
