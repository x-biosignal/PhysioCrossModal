# Multitaper coherence

Computes magnitude-squared coherence using multitaper spectral
estimation with Discrete Prolate Spheroidal Sequences (DPSS / Slepian
tapers). Multitaper estimation provides lower variance than Welch's
method for short signals or when frequency resolution is critical.

## Usage

``` r
multitaperCoherence(
  x,
  y = NULL,
  sr = NULL,
  nw = 4,
  k = NULL,
  freq_range = NULL,
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

- nw:

  Numeric time-bandwidth product (default 4). Controls the trade-off
  between frequency resolution and variance.

- k:

  Integer number of tapers (default `2 * nw - 1`).

- freq_range:

  Optional numeric vector `c(low, high)` to restrict output frequencies.

- modality_x, modality_y:

  Character modality names for MPE input.

- channels_x, channels_y:

  Integer channel indices (default 1).

- ...:

  Currently unused.

## Value

A list with components:

- coherence:

  Numeric vector of magnitude-squared coherence values in \\\[0, 1\]\\.

- frequencies:

  Numeric vector of corresponding frequencies in Hz.

- confidence_limit:

  Numeric scalar giving the 95\\ threshold based on the number of
  tapers.

## References

Thomson, D. J. (1982). Spectrum estimation and harmonic analysis.
*Proceedings of the IEEE*, 70(9), 1055–1096.

Carter, G. C. (1987). Coherence and time delay estimation. *Proceedings
of the IEEE*, 75(2), 236–255.

## See also

[`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md),
[`crossSpectrum()`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossSpectrum.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
sr <- 500
t <- seq(0, 10, length.out = sr * 10)
x <- sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
y <- 0.8 * sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
result <- multitaperCoherence(x, y, sr = sr)
```
