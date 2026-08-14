# Corrected imaginary phase-locking value (ciPLV)

ciPLV (Bruna et al. 2018) corrects the imaginary part of the PLV for the
size of its real part, giving a phase-synchrony estimate that is robust
to zero-lag (volume-conduction) coupling.

## Usage

``` r
ciPLV(
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

A single ciPLV magnitude in `[0, 1]`.

## References

Bruna R, et al. (2018). Phase locking value revisited. *J Neural Eng*
15:056011.
