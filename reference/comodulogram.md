# Comodulogram (phase-frequency by amplitude-frequency PAC)

Computes phase-amplitude coupling over a grid of phase frequencies and
amplitude frequencies, yielding a comodulogram matrix whose peak locates
the dominant coupling.

## Usage

``` r
comodulogram(
  x,
  y = NULL,
  sr = NULL,
  phase_freqs,
  amp_freqs,
  method = c("tort", "canolty", "ozkurt", "plv"),
  phase_bw = 2,
  amp_bw = 10,
  n_bins = 18L,
  order = 4L,
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L
)
```

## Arguments

- x, y, sr, method, n_bins, order, modality_x, modality_y, channels_x,
  channels_y:

  As in
  [`phaseAmplitudeCoupling()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseAmplitudeCoupling.md).

- phase_freqs:

  Numeric vector of phase centre frequencies in Hz.

- amp_freqs:

  Numeric vector of amplitude centre frequencies in Hz.

- phase_bw:

  Half-bandwidth of each phase band in Hz (default: 2).

- amp_bw:

  Half-bandwidth of each amplitude band in Hz (default: 10).

## Value

A list with the `matrix` (rows = phase frequencies, columns = amplitude
frequencies), `phase_freqs`, `amp_freqs`, `method`, and `peak` (the
(phase, amplitude) frequency of the maximum and its value).

## See also

[`phaseAmplitudeCoupling()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseAmplitudeCoupling.md),
[`plotComodulogram()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotComodulogram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1); t <- (0:4999) / 500; ph <- 2 * pi * 6 * t
sig <- sin(ph) + (0.1 + (1 + cos(ph - pi)) / 2) * sin(2 * pi * 50 * t)
cm <- comodulogram(sig, sr = 500, phase_freqs = c(4, 6, 8),
                   amp_freqs = c(30, 50, 70))
cm$peak
} # }
```
