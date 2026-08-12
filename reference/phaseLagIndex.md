# Phase Lag Index (PLI)

Computes the Phase Lag Index between two signals. PLI measures the
asymmetry of the distribution of phase differences and is insensitive to
zero-lag (volume conduction) effects.

## Usage

``` r
phaseLagIndex(
  x,
  y = NULL,
  sr = NULL,
  freq_band,
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

- pli:

  Numeric scalar in \[0, 1\] – the Phase Lag Index.

- phase_diff:

  Numeric vector of instantaneous phase differences (radians).

## Details

\$\$\text{PLI} = \left\|\frac{1}{N} \sum\_{t=1}^{N}
\text{sign}(\sin(\phi_x(t) - \phi_y(t)))\right\|\$\$

A PLI of 0 indicates either no coupling or symmetric (zero-lag)
coupling. A PLI of 1 indicates perfect non-zero-lag phase locking.

## References

Stam, C. J., Nolte, G., & Daffertshofer, A. (2007). Phase lag index:
Assessment of functional connectivity from multi channel EEG and MEG
with diminished bias from common sources. *Human Brain Mapping*, 28(11),
1178–1193.

## See also

[`phaseLockingValue()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLockingValue.md),
[`weightedPLI()`](https://x-biosignal.github.io/PhysioCrossModal/reference/weightedPLI.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
sr <- 500
t <- seq(0, 2, length.out = sr * 2)
x <- sin(2 * pi * 10 * t)
y <- sin(2 * pi * 10 * t + pi / 3)
result <- phaseLagIndex(x, y, sr = sr, freq_band = c(8, 12))
result$pli
#> [1] 0.99
```
