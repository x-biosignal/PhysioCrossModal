# Elastic (SRVF) alignment of two curves

Aligns curve `y` to curve `x` by an optimal time warping computed in the
Square-Root Velocity Function representation, separating phase (the
warping) from amplitude (the residual shape difference). A pure-R
dynamic-programming solver is used by default; if the fdasrvf package is
available it is used for higher accuracy.

## Usage

``` r
elasticAlign(
  x,
  y,
  t = NULL,
  lambda = 0,
  method = c("srvf"),
  max_step = 6L,
  smooth = FALSE,
  use_fdasrvf = NA
)
```

## Arguments

- x, y:

  Numeric curves of the same length.

- t:

  Optional grid (defaults to an evenly spaced grid on 0..1).

- lambda:

  Non-negative warping-roughness penalty (passed to fdasrvf when used;
  the pure-R solver ignores it beyond boundary/monotonicity).

- method:

  Currently `"srvf"`.

- max_step:

  Maximum DP step (bounds the local warping slope) for the pure-R solver
  (default: 6).

- smooth:

  Smooth the pure-R warping with a monotone spline (default: `FALSE`).

- use_fdasrvf:

  `NA` (default) uses fdasrvf if installed; `TRUE` requires it; `FALSE`
  forces the pure-R solver.

## Value

A list with `gamma` (the warping), `aligned` (`y` warped to `x`),
`amplitude_distance`, `phase_distance`, the SRVFs, and the solver used.

## References

Srivastava, A., Wu, W., Kurtek, S., Klassen, E., & Marron, J. S. (2011).
Registration of functional data using Fisher-Rao metric.
arXiv:1103.3817.

## See also

[`srvfMean()`](https://x-biosignal.github.io/PhysioCrossModal/reference/srvfMean.md),
[`warpApply()`](https://x-biosignal.github.io/PhysioCrossModal/reference/warpApply.md)

## Examples

``` r
t <- seq(0, 1, length.out = 80)
x <- exp(-((t - 0.5)^2) / (2 * 0.08^2))
y <- exp(-((t - 0.65)^2) / (2 * 0.08^2))
al <- elasticAlign(x, y)
al$amplitude_distance
#> [1] 0.06749576
```
