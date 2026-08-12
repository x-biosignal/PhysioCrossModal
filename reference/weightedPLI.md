# Weighted Phase Lag Index (wPLI)

Computes the weighted Phase Lag Index between two signals. wPLI weights
each phase-difference sample by the magnitude of the imaginary component
of the cross-spectrum, reducing the influence of noise sources that
produce phase differences near 0 or pi.

## Usage

``` r
weightedPLI(
  x,
  y = NULL,
  sr = NULL,
  freq_band,
  debiased = TRUE,
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

- debiased:

  Logical; if TRUE (default), also compute the debiased wPLI estimator.

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

- wpli:

  Numeric scalar in \[0, 1\] – the weighted PLI.

- wpli_debiased:

  Numeric scalar – the debiased wPLI, or NULL if `debiased = FALSE`.

- n_samples:

  Integer – number of time-domain samples used.

## Details

\$\$\text{wPLI} = \frac{\left\|\sum \|\text{Im}(S\_{xy})\| \cdot
\text{sign}(\text{Im}(S\_{xy}))\right\|}{\sum
\|\text{Im}(S\_{xy})\|}\$\$

When `debiased = TRUE`, the debiased estimator from Vinck et al. (2011)
is also returned: \$\$\text{wPLI}\_{\text{debiased}} = \frac{(\sum
\text{Im}(S\_{xy}))^2 - \sum \text{Im}(S\_{xy})^2}{(\sum
\|\text{Im}(S\_{xy})\|)^2 - \sum \text{Im}(S\_{xy})^2}\$\$

## References

Vinck, M., Oostenveld, R., van Wingerden, M., Battaglia, F., & Pennartz,
C. M. A. (2011). An improved index of phase-synchronization for
electrophysiological data in the presence of volume-conduction, noise
and sample-size bias. *NeuroImage*, 55(4), 1548–1565.

Lachaux, J.-P., Rodriguez, E., Martinerie, J., & Varela, F. J. (1999).
Measuring phase synchrony in brain signals. *Human Brain Mapping*, 8(4),
194–208.

## See also

[`phaseLockingValue()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLockingValue.md),
[`phaseLagIndex()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLagIndex.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
sr <- 500
t <- seq(0, 2, length.out = sr * 2)
x <- sin(2 * pi * 10 * t)
y <- sin(2 * pi * 10 * t + pi / 4)
result <- weightedPLI(x, y, sr = sr, freq_band = c(8, 12))
result$wpli
#> [1] 0.999989
```
