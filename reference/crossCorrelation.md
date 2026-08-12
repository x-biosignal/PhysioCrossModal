# Cross-correlation between two signals

Computes the cross-correlation between two signals at various lags,
measuring their time-domain coupling. Cross-correlation quantifies how
similar one signal is to a time-shifted version of another. A positive
peak lag indicates that `y` leads `x` (i.e., `y` must be shifted forward
to align with `x`), while a negative peak lag indicates `x` leads `y`.

## Usage

``` r
crossCorrelation(
  x,
  y = NULL,
  sr = NULL,
  max_lag = NULL,
  normalize = TRUE,
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

- max_lag:

  Integer maximum lag in samples. If NULL, defaults to
  `floor(length(x) / 4)`.

- normalize:

  Logical; if TRUE (default), compute normalized cross-correlation
  (Pearson-like, values in \[-1, 1\]).

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

- correlation:

  Numeric vector of cross-correlation values at each lag.

- lags:

  Integer vector of lag values in samples.

- lag_seconds:

  Numeric vector of lag values in seconds.

- peak_lag:

  Integer lag (in samples) at which the absolute correlation is
  maximised.

- peak_lag_seconds:

  Numeric peak lag converted to seconds.

- peak_correlation:

  Numeric cross-correlation value at the peak lag.

## Details

The function accepts numeric vectors, PhysioExperiment objects, or a
MultiPhysioExperiment with named modalities, using
`.extract_signal_pair()` internally for flexible input handling.

## References

Chatfield, C. (2004). *The Analysis of Time Series: An Introduction*
(6th ed.). Chapman & Hall/CRC.

## See also

[`slidingCrossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/slidingCrossCorrelation.md),
[`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
# Simple example: y is a delayed copy of x
sr <- 500
set.seed(1)
x <- rnorm(1000)
y <- c(rep(0, 10), x[1:990])
result <- crossCorrelation(x, y, sr = sr)
result$peak_lag       # should be 10
#> [1] -10
result$peak_lag_seconds
#> [1] -0.02
```
