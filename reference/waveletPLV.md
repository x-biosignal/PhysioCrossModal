# Time-frequency wavelet PLV

Computes Phase Locking Value in the time-frequency domain using complex
Morlet wavelets. The phase difference between the two signals is
computed at each time-frequency point, and PLV is the magnitude of the
smoothed unit-phase vector:

## Usage

``` r
waveletPLV(
  x,
  y = NULL,
  sr = NULL,
  frequencies = seq(1, 40, by = 1),
  n_cycles = 7,
  smoothing_cycles = 3,
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L,
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

- frequencies:

  Numeric vector of centre frequencies in Hz (default
  `seq(1, 40, by = 1)`).

- n_cycles:

  Numeric number of wavelet cycles (default 7).

- smoothing_cycles:

  Numeric number of cycles for the temporal smoothing Gaussian (default
  3).

- modality_x, modality_y:

  Character modality names for MPE input.

- channels_x, channels_y:

  Integer channel indices (default 1).

- ...:

  Currently unused.

## Value

A list with components:

- plv:

  Numeric matrix `[time x frequency]` of PLV values in \\\[0, 1\]\\.

- frequencies:

  Numeric vector of centre frequencies.

- times:

  Numeric vector of time points (seconds from start).

- coi:

  Numeric vector of Cone of Influence frequencies. At each time point,
  frequencies below this value are affected by edge artifacts.

## Details

\$\$\text{PLV}(t,f) = \left\|\langle e^{i\Delta\phi(t,f)}
\rangle\right\|\$\$

## References

Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet
analysis. *Bulletin of the American Meteorological Society*, 79(1),
61–78.

Lachaux, J.-P., Rodriguez, E., Martinerie, J., & Varela, F. J. (1999).
Measuring phase synchrony in brain signals. *Human Brain Mapping*, 8(4),
194–208.

## See also

[`waveletCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletCoherence.md),
[`phaseLockingValue()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLockingValue.md),
[`plotWaveletCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotWaveletCoherence.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
sr <- 200
t <- seq(0, 2, length.out = sr * 2)
x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
y <- sin(2 * pi * 10 * t + pi/4) + 0.3 * rnorm(length(t))
result <- waveletPLV(x, y, sr = sr, frequencies = seq(5, 20))
```
