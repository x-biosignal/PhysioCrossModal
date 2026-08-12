# Sliding-window cross-correlation

Computes cross-correlation in sliding (overlapping) windows to track how
time-domain coupling varies over time. For each window position,
[`crossCorrelation`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossCorrelation.md)
is called and the results are assembled into a matrix.

## Usage

``` r
slidingCrossCorrelation(
  x,
  y = NULL,
  sr = NULL,
  window_sec = 1,
  step_sec = 0.5,
  max_lag = NULL,
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L
)
```

## Arguments

- x:

  Numeric vector, PhysioExperiment, or MultiPhysioExperiment.

- y:

  Numeric vector or PhysioExperiment, or NULL when `x` is a
  MultiPhysioExperiment.

- sr:

  Numeric sampling rate in Hz (required when x/y are numeric).

- window_sec:

  Numeric window length in seconds (default 1).

- step_sec:

  Numeric step size in seconds (default 0.5).

- max_lag:

  Integer maximum lag in samples for each window. If NULL, defaults to
  `floor(window_samples / 4)`.

- modality_x:

  Character modality name in MPE for the x signal.

- modality_y:

  Character modality name in MPE for the y signal.

- channels_x:

  Integer which channel to extract from x (default 1).

- channels_y:

  Integer which channel to extract from y (default 1).

## Value

A named list with components:

- correlations:

  Numeric matrix of dimensions (n_windows x n_lags) containing
  cross-correlation values.

- times:

  Numeric vector of window centre times in seconds.

- lags:

  Integer vector of lag values in samples.

- peak_lags:

  Numeric vector of peak lags (in samples) for each window.

- peak_correlations:

  Numeric vector of peak correlation values for each window.

## References

Chatfield, C. (2004). *The Analysis of Time Series: An Introduction*
(6th ed.). Chapman & Hall/CRC.

## See also

[`crossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossCorrelation.md),
[`plotCouplingTimecourse()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotCouplingTimecourse.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
sr <- 500
set.seed(1)
x <- rnorm(5000)
y <- c(rep(0, 10), x[1:4990])
result <- slidingCrossCorrelation(x, y, sr = sr,
                                   window_sec = 1, step_sec = 0.5)
dim(result$correlations)
#> [1]  19 251
result$peak_lags
#>  [1] -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10
```
