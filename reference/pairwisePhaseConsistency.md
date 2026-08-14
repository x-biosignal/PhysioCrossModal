# Pairwise phase consistency (PPC)

The bias-free phase-synchrony estimator of Vinck et al. (2010): \\PPC =
(\|\sum e^{i\Delta\phi}\|^2 - N) / (N(N-1))\\. Unlike the PLV it is not
inflated by the sample size.

## Usage

``` r
pairwisePhaseConsistency(
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

A single PPC value (approaches 1 for perfect locking, 0 for none).

## References

Vinck M, et al. (2010). The pairwise phase consistency. *NeuroImage*
51:112-122.
