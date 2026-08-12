# Phase Locking Value (PLV)

Computes the Phase Locking Value between two signals, measuring the
consistency of phase difference across time. A PLV of 1 indicates
perfect phase locking; a PLV of 0 indicates no consistent phase
relationship.

## Usage

``` r
phaseLockingValue(
  x,
  y = NULL,
  sr = NULL,
  freq_band,
  method = c("hilbert", "wavelet"),
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

  Numeric vector or PhysioExperiment (NULL when `x` is an MPE).

- sr:

  Numeric sampling rate in Hz (required when x/y are numeric).

- freq_band:

  Numeric vector of length 2 specifying the frequency band
  `c(low, high)` in Hz. Required.

- method:

  Character; phase extraction method. Currently only `"hilbert"` is
  supported. `"wavelet"` is reserved for future use.

- modality_x:

  Character modality name in MPE for x signal.

- modality_y:

  Character modality name in MPE for y signal.

- channels_x:

  Integer channel index to extract from x (default 1).

- channels_y:

  Integer channel index to extract from y (default 1).

## Value

A named list with components:

- plv:

  Numeric scalar in \[0, 1\] – the Phase Locking Value.

- phase_diff:

  Numeric vector of instantaneous phase differences (radians, \[-pi,
  pi\]).

## Details

The signals are bandpass-filtered to `freq_band` and then Hilbert-
transformed to extract instantaneous phase. PLV is computed as:
\$\$\text{PLV} = \left\|\frac{1}{N} \sum\_{t=1}^{N} e^{i(\phi_x(t) -
\phi_y(t))}\right\|\$\$

## References

Lachaux, J.-P., Rodriguez, E., Martinerie, J., & Varela, F. J. (1999).
Measuring phase synchrony in brain signals. *Human Brain Mapping*, 8(4),
194–208.

## See also

[`phaseLagIndex()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLagIndex.md),
[`weightedPLI()`](https://x-biosignal.github.io/PhysioCrossModal/reference/weightedPLI.md),
[`waveletPLV()`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletPLV.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
# Two phase-locked sinusoids
sr <- 500
t <- seq(0, 2, length.out = sr * 2)
x <- sin(2 * pi * 10 * t)
y <- sin(2 * pi * 10 * t + pi / 4)
result <- phaseLockingValue(x, y, sr = sr, freq_band = c(8, 12))
result$plv
#> [1] 0.993784
```
