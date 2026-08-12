# Lagged coherence between two signals

Convenience wrapper for
[`coherence`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md)
with `type = "lagged"`, the lagged (non-instantaneous) coherence of
Pascual-Marqui (2007), `Im(C) / sqrt(1 - Re(C)^2)` where `C` is the
complex coherency. Like the imaginary part of coherency it is
insensitive to zero-lag / volume-conduction coupling.

## Usage

``` r
laggedCoherence(
  x,
  y = NULL,
  sr = NULL,
  freq_range = NULL,
  nperseg = 256L,
  noverlap = NULL,
  method = c("welch", "multitaper"),
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

  Numeric vector or PhysioExperiment, or NULL when `x` is an MPE.

- sr:

  Numeric sampling rate in Hz (required when x/y are numeric).

- freq_range:

  Numeric vector of length 2 giving the frequency range `c(low, high)`
  in Hz to restrict the output. NULL returns all frequencies.

- nperseg:

  Integer segment length for Welch's method (default 256).

- noverlap:

  Integer overlap between segments, or NULL for `floor(nperseg / 2)`.

- method:

  Character spectral estimation method. Currently only `"welch"` is
  implemented; `"multitaper"` is reserved for future use.

- modality_x:

  Character modality name in MPE for the x signal.

- modality_y:

  Character modality name in MPE for the y signal.

- channels_x:

  Integer which channel to extract from x (default 1).

- channels_y:

  Integer which channel to extract from y (default 1).

## Value

A named list as returned by
[`coherence`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md),
with `type = "lagged"`.

## References

Pascual-Marqui, R. D. (2007). Instantaneous and lagged measurements of
linear and nonlinear dependence between groups of multivariate time
series. *arXiv:0711.1455*.

## See also

[`coherence`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md)

## Examples

``` r
sr <- 200
t <- seq(0, 10, length.out = sr * 10)
x <- sin(2 * pi * 10 * t)
y <- sin(2 * pi * 10 * t + pi / 2)          # 90-degree lag
lc <- laggedCoherence(x, y, sr = sr, freq_range = c(8, 12))
max(abs(lc$coherence))
#> [1] 1
```
