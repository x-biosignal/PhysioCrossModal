# Phase-amplitude coupling

Quantifies how the amplitude of a fast oscillation (in `amp_band`) is
modulated by the phase of a slow oscillation (in `phase_band`). The
phase is taken from `x` and the amplitude envelope from `y` (or from `x`
itself when `y` is `NULL`), each via a bandpass filter and the Hilbert
transform. Supported measures are the Tort modulation index (a
normalized Kullback-Leibler distance of the phase-amplitude distribution
from uniform, in 0..1), the Canolty mean vector length (raw), the Ozkurt
normalized mean vector length (in 0..1), and a phase-locking value.

## Usage

``` r
phaseAmplitudeCoupling(
  x,
  y = NULL,
  sr = NULL,
  phase_band = c(4, 8),
  amp_band = c(30, 80),
  method = c("tort", "canolty", "ozkurt", "plv"),
  n_bins = 18L,
  order = 4L,
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L
)
```

## Arguments

- x:

  A numeric signal or PhysioExperiment (phase source), or a
  MultiPhysioExperiment.

- y:

  A numeric signal or PhysioExperiment (amplitude source), or `NULL` to
  use `x` for both.

- sr:

  Sampling rate in Hz (required for numeric input).

- phase_band:

  Numeric length-2 phase frequency band in Hz (default: `c(4, 8)`).

- amp_band:

  Numeric length-2 amplitude frequency band in Hz (default:
  `c(30, 80)`).

- method:

  `"tort"`, `"canolty"`, `"ozkurt"`, or `"plv"`.

- n_bins:

  Number of phase bins for the Tort modulation index (default: 18).

- order:

  Bandpass filter order (default: 4).

- modality_x, modality_y, channels_x, channels_y:

  Passed to the shared signal-pair extraction.

## Value

A list with `pac` (the coupling value), `method`, `phase_band`,
`amp_band`, and, for the Tort method, the per-bin normalized amplitude
`distribution` and `bin_centers`.

## References

Tort, A. B. L., et al. (2010). Measuring phase-amplitude coupling
between neuronal oscillations of different frequencies. Journal of
Neurophysiology, 104(2), 1195-1210.

Canolty, R. T., et al. (2006). High gamma power is phase-locked to theta
oscillations in human neocortex. Science, 313(5793), 1626-1628.

## See also

[`comodulogram()`](https://x-biosignal.github.io/PhysioCrossModal/reference/comodulogram.md),
[`plotComodulogram()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotComodulogram.md),
[`surrogateTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateTest.md)

## Examples

``` r
set.seed(1)
t <- (0:4999) / 500
ph <- 2 * pi * 6 * t
sig <- sin(ph) + (0.1 + (1 + cos(ph - pi)) / 2) * sin(2 * pi * 50 * t)
phaseAmplitudeCoupling(sig, sr = 500, phase_band = c(4, 8),
                       amp_band = c(35, 65))$pac
#> [1] 0.06449726
```
