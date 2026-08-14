# Partial Directed Coherence for multichannel / multimodal signals

Frequency-resolved directed connectivity via Partial Directed Coherence
(Baccala & Sameshima 2001), computed from the coefficient matrix of an
MVAR model
([`PhysioCore::mvarFit()`](https://x-biosignal.github.io/PhysioCore//reference/mvarFit.html));
numerically identical to `PhysioEEG`'s `eegPDC`. Unlike DTF, PDC
reflects only direct influences, so a purely indirect pathway gives PDC
near zero.

## Usage

``` r
partialDirectedCoherence(
  x,
  y = NULL,
  sr = NULL,
  order = 5L,
  freqs = NULL,
  generalized = TRUE,
  band = NULL,
  mvar_method = "ols",
  modality_x = NULL,
  modality_y = NULL,
  channels_x = NULL,
  channels_y = NULL
)
```

## Arguments

- x:

  A multichannel matrix (time x channels) or PhysioExperiment (with
  `y = NULL`), or one signal / PhysioExperiment of a pair.

- y:

  The second signal / PhysioExperiment of a pair, or `NULL` for a single
  multichannel input.

- sr:

  Sampling rate in Hz (required for numeric-matrix input).

- order:

  MVAR model order (default: 5), or `NULL` to select automatically.

- freqs:

  Frequencies in Hz (default: 128 points from 0 to the Nyquist
  frequency).

- generalized:

  Use generalized PDC (inverse-noise weighting; default `TRUE`).

- band:

  Optional length-2 band in Hz over which to average the returned matrix
  (default: all frequencies).

- mvar_method:

  MVAR estimator (default: `"ols"`).

- modality_x, modality_y, channels_x, channels_y:

  Passed to the shared signal-pair extraction when `y` is supplied.

## Value

A list as in
[`directedTransferFunction()`](https://x-biosignal.github.io/PhysioCrossModal/reference/directedTransferFunction.md)
with a `pdc` array and band-averaged `matrix`; for a two-channel input,
`pdc_xy` and `pdc_yx`.

## References

Baccala, L. A., & Sameshima, K. (2001). Partial directed coherence: a
new concept in neural structure determination. Biological Cybernetics,
84(6), 463-474.

## See also

[`directedTransferFunction()`](https://x-biosignal.github.io/PhysioCrossModal/reference/directedTransferFunction.md),
[`transferEntropy()`](https://x-biosignal.github.io/PhysioCrossModal/reference/transferEntropy.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(1500), ncol = 3)
res <- partialDirectedCoherence(X, sr = 100, order = 4)
dim(res$pdc)
#> [1]   3   3 128
```
