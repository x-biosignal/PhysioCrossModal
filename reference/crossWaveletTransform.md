# Cross-wavelet transform (XWT) with red-noise significance

Computes the cross-wavelet transform \\W\_{xy} = W_x W_y^\*\\ of two
signals via the Morlet wavelet: its modulus is the common time-frequency
power and its argument is the phase lead/lag between the signals. Common
power is tested against an AR(1) red-noise background (Torrence & Compo
1998; Grinsted et al. 2004) by Monte-Carlo: AR(1) surrogates matched to
each signal's lag-1 autocorrelation and variance are transformed
identically, and the per- frequency `siglvl` quantile of the surrogate
cross-power sets the threshold.

## Usage

``` r
crossWaveletTransform(
  x,
  y = NULL,
  sr = NULL,
  frequencies = seq(1, 40, by = 1),
  n_cycles = 7,
  significance = TRUE,
  n_surrogates = 100L,
  siglvl = 0.95,
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

- frequencies:

  Frequencies (Hz) to analyse (default `seq(1, 40, 1)`).

- n_cycles:

  Morlet wavelet cycles (default 7).

- significance:

  If `TRUE` (default), compute red-noise significance.

- n_surrogates:

  Number of AR(1) surrogate pairs (default 100).

- siglvl:

  Significance quantile (default 0.95).

- modality_x, modality_y, channels_x, channels_y:

  Passed to the pair extractor.

## Value

An object of class `"cross_wavelet"`: a list with `power`
(`time x frequency` cross-wavelet modulus), `phase` (relative phase,
radians; positive = `x` leads `y`), `cross` (complex XWT),
`frequencies`, `times`, `coi` (cone of influence), and — when
`significance` — `sig_level` (per frequency), `sig_ratio` (power /
level), and `significant` (logical mask).

## References

Torrence C, Compo GP (1998). A practical guide to wavelet analysis.
*Bull Amer Meteor Soc* 79:61-78. Grinsted A, Moore JC, Jevrejeva S
(2004). Application of the cross wavelet transform and wavelet coherence
to geophysical time series. *Nonlin Processes Geophys* 11:561-566.

## Examples

``` r
set.seed(1); sr <- 100; t <- seq(0, 20, by = 1 / sr)
x <- sin(2 * pi * 10 * t) + rnorm(length(t), sd = 0.5)
y <- sin(2 * pi * 10 * t - pi / 2) + rnorm(length(t), sd = 0.5)  # y lags x
xw <- crossWaveletTransform(x, y, sr = sr, frequencies = seq(4, 16, by = 1),
                            n_surrogates = 50)
mean(xw$significant)
#> [1] 0.264906
```
