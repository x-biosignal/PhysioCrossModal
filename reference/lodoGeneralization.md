# Leave-one-site-out generalization benchmark

Evaluates model transportability by training on all-but-one sites and
testing on the held-out site (LODO). Supports continuous outcomes
(`family = "gaussian"`) and binary outcomes (`family = "binomial"`).

## Usage

``` r
lodoGeneralization(
  data,
  outcome,
  site,
  features = NULL,
  family = c("gaussian", "binomial"),
  positive_class = NULL,
  threshold = 0.5,
  min_train_rows = 20L,
  scale_features = TRUE
)
```

## Arguments

- data:

  Data frame containing outcome, site label, and features.

- outcome:

  Outcome column name.

- site:

  Site/facility column name used for LODO splitting.

- features:

  Feature column names. If `NULL`, numeric columns excluding `outcome`
  and `site` are used.

- family:

  Modeling family: `"gaussian"` or `"binomial"`.

- positive_class:

  Positive class label for binomial metrics. If `NULL`, the second
  factor level is used.

- threshold:

  Classification threshold for binomial predictions.

- min_train_rows:

  Minimum training rows required per fold.

- scale_features:

  Logical; z-score features using training statistics.

## Value

A list with:

- fold_metrics:

  Per-site metrics.

- predictions:

  Row-level held-out predictions.

- aggregate:

  Mean metrics across folds.

- settings:

  Benchmark settings and feature list.

## See also

[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md),
[`surrogateTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateTest.md)

## Examples

``` r
set.seed(1)
n <- 120
df <- data.frame(
  site_id = rep(c("A", "B", "C"), each = 40),
  f1 = rnorm(n),
  f2 = rnorm(n)
)
df$y <- 0.8 * df$f1 - 0.4 * df$f2 + rnorm(n, sd = 0.2)

res <- lodoGeneralization(
  data = df,
  outcome = "y",
  site = "site_id",
  family = "gaussian"
)
head(res$fold_metrics)
#>   holdout_site n_train n_test status      rmse       mae        r2
#> 1            A      80     40     ok 0.2054895 0.1611881 0.9353797
#> 2            B      80     40     ok 0.1707508 0.1387755 0.9730365
#> 3            C      80     40     ok 0.2024607 0.1644832 0.9270859
```
