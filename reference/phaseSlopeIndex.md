# Phase Slope Index (directed, volume-conduction robust)

The Phase Slope Index (Nolte et al. 2008) estimates the direction of
information flow from the slope of the coherency phase across a
frequency band, and is insensitive to instantaneous (zero-lag) mixing.
Positive values indicate that `x` leads (drives) `y`.

## Usage

``` r
phaseSlopeIndex(
  x,
  y = NULL,
  sr = NULL,
  freq_band,
  nperseg = 256L,
  noverlap = NULL,
  normalize = TRUE,
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L
)
```

## Arguments

- x, y:

  Signals (or a single object from which a pair is extracted).

- sr:

  Sampling rate in Hz.

- freq_band:

  Numeric `c(low, high)` band in Hz.

- nperseg, noverlap:

  Welch segment length / overlap for the cross-spectrum.

- normalize:

  If `TRUE` (default), scale into `[-1, 1]` by the summed coherency
  magnitude; if `FALSE`, return the raw imaginary sum.

- modality_x, modality_y, channels_x, channels_y:

  Passed to the internal signal-pair extractor.

## Value

A single Phase Slope Index value.

## References

Nolte G, et al. (2008). Robustly estimating the flow direction of
information in complex physical systems. *Phys Rev Lett* 100:234101.
[doi:10.1103/PhysRevLett.100.234101](https://doi.org/10.1103/PhysRevLett.100.234101)

## Examples

``` r
set.seed(1); sr <- 200; n <- 4000
x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2)); x[is.na(x)] <- 0
y <- c(rep(0, 8), x[seq_len(n - 8)])          # y lags x
phaseSlopeIndex(x, y, sr = sr, freq_band = c(1, 40))
#> [1] 0.1940465
```
