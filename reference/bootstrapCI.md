# Bootstrap confidence interval for coupling

Computes a bootstrap confidence interval for the coupling statistic
between two signals using the moving-block bootstrap (to preserve
temporal autocorrelation).

## Usage

``` r
bootstrapCI(
  x,
  y = NULL,
  sr = NULL,
  method,
  n_boot = 199L,
  ci = 0.95,
  block_len = NULL,
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L,
  cores = 1L,
  ...
)
```

## Arguments

- x:

  Numeric vector, PhysioExperiment, or MultiPhysioExperiment.

- y:

  Numeric vector or PhysioExperiment, or NULL when `x` is an MPE.

- sr:

  Numeric sampling rate in Hz (required when x/y are numeric).

- method:

  Character coupling method (same options as
  [`surrogateTest`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateTest.md)).

- n_boot:

  Integer number of bootstrap replicates (default 199).

- ci:

  Numeric confidence level (default 0.95).

- block_len:

  Integer block length for block bootstrap, or NULL for automatic
  (`ceiling(sqrt(n))`).

- modality_x, modality_y:

  Character modality names for MPE input.

- channels_x, channels_y:

  Integer channel indices (default 1).

- cores:

  Integer number of parallel cores to use (default `1L`). When
  `cores > 1`, bootstrap replicates are computed in parallel.

- ...:

  Additional arguments passed to the coupling function.

## Value

A list with components:

- observed:

  The full coupling result from the original signals.

- statistic:

  Numeric scalar: the extracted coupling statistic.

- ci_lower:

  Numeric scalar: lower CI bound.

- ci_upper:

  Numeric scalar: upper CI bound.

- ci_level:

  Numeric scalar: confidence level used.

- boot_distribution:

  Numeric vector of bootstrap statistics.

## References

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall/CRC.

## See also

[`surrogateTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateTest.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
sr <- 500
t <- seq(0, 2, length.out = sr * 2)
x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
y <- 0.8 * sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
result <- bootstrapCI(x, y, sr = sr, method = "coherence",
                      n_boot = 19, nperseg = 128L)
c(result$ci_lower, result$ci_upper)
#> [1] 0.9705681 0.9873852
```
