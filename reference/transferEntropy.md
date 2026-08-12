# Transfer entropy (model-free directed information)

Estimates the transfer entropy from `x` to `y` and back (Schreiber
2000), a model-free measure of directed information flow, with a
histogram (binning) estimator or the Kraskov-Stoegbauer-Grassberger
(KSG) k-nearest-neighbour estimator. With `effective = TRUE` the
estimate is bias-corrected by subtracting the mean transfer entropy of
IAAFT surrogates of the source, so independent signals give an effective
transfer entropy near zero.

## Usage

``` r
transferEntropy(
  x,
  y = NULL,
  sr = NULL,
  k = 1L,
  l = 1L,
  delay = 1L,
  estimator = c("ksg", "histogram"),
  knn = 4L,
  bins = 8L,
  effective = FALSE,
  n_surrogate = 19L,
  seed = NULL,
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L
)
```

## Arguments

- x, y:

  Two numeric signals or PhysioExperiment objects (or a
  MultiPhysioExperiment via `x`).

- sr:

  Sampling rate in Hz (required for numeric input).

- k:

  Target (destination) embedding length (default: 1).

- l:

  Source embedding length (default: 1).

- delay:

  Source-to-target lag in samples (default: 1).

- estimator:

  `"ksg"` (k-NN, in nats) or `"histogram"` (binning, in bits).

- knn:

  Number of neighbours for the KSG estimator (default: 4).

- bins:

  Number of bins per dimension for the histogram estimator (default: 8).

- effective:

  Subtract the surrogate (IAAFT) null mean (default: `FALSE`).

- n_surrogate:

  Number of surrogates for `effective` (default: 19).

- seed:

  Optional RNG seed for reproducible surrogates.

- modality_x, modality_y, channels_x, channels_y:

  Passed to the shared signal-pair extraction.

## Value

A list with `te_xy`, `te_yx`, and `net` (`te_xy - te_yx`); when
`effective = TRUE` also `eff_xy`, `eff_yx`, and the surrogate null
standard deviations.

## References

Schreiber, T. (2000). Measuring information transfer. Physical Review
Letters, 85(2), 461-464.

Kraskov, A., Stoegbauer, H., & Grassberger, P. (2004). Estimating mutual
information. Physical Review E, 69(6), 066138.

## See also

[`directedTransferFunction()`](https://x-biosignal.github.io/PhysioCrossModal/reference/directedTransferFunction.md),
[`grangerCausality()`](https://x-biosignal.github.io/PhysioCrossModal/reference/grangerCausality.md),
[`surrogateTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateTest.md)

## Examples

``` r
set.seed(1)
x <- rnorm(600); y <- c(rep(0, 5), 0.7 * x[1:595]) + rnorm(600, sd = 0.5)
transferEntropy(x, y, sr = 100, delay = 5, estimator = "histogram")$net
#> [1] 0.5250553
```
