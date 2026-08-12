# Time-frequency wavelet coherence

Computes wavelet coherence between two signals using complex Morlet
wavelets. The cross-wavelet spectrum and auto-spectra are smoothed with
a Gaussian temporal window (width proportional to
`smoothing_cycles / frequency`), and coherence is computed as:

## Usage

``` r
waveletCoherence(
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

- coherence:

  Numeric matrix `[time x frequency]` of coherence values in \\\[0,
  1\]\\.

- phase:

  Numeric matrix `[time x frequency]` of phase differences (radians).

- frequencies:

  Numeric vector of centre frequencies.

- times:

  Numeric vector of time points (seconds from start).

- coi:

  Numeric vector of Cone of Influence frequencies. At each time point,
  frequencies below this value are affected by edge artifacts.

## Details

\$\$C(t,f) = \frac{\|\langle W\_{xy}(t,f) \rangle\|^2}{\langle
\|W_x(t,f)\|^2 \rangle \cdot \langle \|W_y(t,f)\|^2 \rangle}\$\$

## References

Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet
analysis. *Bulletin of the American Meteorological Society*, 79(1),
61–78.

Grinsted, A., Moore, J. C., & Jevrejeva, S. (2004). Application of the
cross wavelet transform and wavelet coherence to geophysical time
series. *Nonlinear Processes in Geophysics*, 11(5/6), 561–566.

## See also

[`waveletPLV()`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletPLV.md),
[`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md),
[`plotWaveletCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotWaveletCoherence.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
sr <- 200
t <- seq(0, 2, length.out = sr * 2)
x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
y <- 0.8 * sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
result <- waveletCoherence(x, y, sr = sr, frequencies = seq(5, 20))
```
