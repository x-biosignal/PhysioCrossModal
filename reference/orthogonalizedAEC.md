# Orthogonalised amplitude-envelope correlation (leakage-corrected)

Amplitude-envelope correlation after pairwise orthogonalisation (Hipp et
al. 2012; Brookes et al. 2012): the instantaneous zero-lag (leakage)
component is projected out before correlating band-limited amplitude
envelopes, so shared instantaneous mixing does not create spurious
coupling. Symmetrised over the two orthogonalisation directions.

## Usage

``` r
orthogonalizedAEC(
  x,
  y = NULL,
  sr = NULL,
  freq_band,
  order = 4L,
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

- freq_band:

  Numeric `c(low, high)` band in Hz.

- order:

  Butterworth filter order for the band-pass (default 4).

- modality_x, modality_y, channels_x, channels_y:

  Passed to the internal signal-pair extractor.

## Value

A single leakage-corrected envelope correlation in `[-1, 1]`.

## References

Hipp JF, et al. (2012). Large-scale cortical correlation structure of
spontaneous oscillatory activity. *Nat Neurosci* 15:884-890.

## Examples

``` r
set.seed(1); sr <- 200; n <- 4000
x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2)); x[is.na(x)] <- 0
y <- c(rep(0, 10), x[seq_len(n - 10)]) + rnorm(n, sd = 0.5)
orthogonalizedAEC(x, y, sr = sr, freq_band = c(8, 12))
#> [1] 0.1046873
```
