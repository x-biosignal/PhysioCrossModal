# Directed Transfer Function for multichannel / multimodal signals

Frequency-resolved directed connectivity via the Directed Transfer
Function (Kaminski & Blinowska 1991), computed from the transfer
function of a multivariate autoregressive model fitted with
[`PhysioCore::mvarFit()`](https://x-biosignal.r-universe.dev/PhysioCore/reference/mvarFit.html).
This cross-modal implementation is numerically identical to
`PhysioEEG`'s `eegDTF` on the same input. Because the transfer function
inverts the whole VAR system, DTF reflects both direct and indirect
pathways.

## Usage

``` r
directedTransferFunction(
  x,
  y = NULL,
  sr = NULL,
  order = 5L,
  freqs = NULL,
  normalized = TRUE,
  ffDTF = FALSE,
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

- normalized:

  Row-normalize the DTF (default: `TRUE`).

- ffDTF:

  Use the full-frequency DTF normalization (default: `FALSE`).

- band:

  Optional length-2 band in Hz over which to average the returned matrix
  (default: all frequencies).

- mvar_method:

  MVAR estimator (default: `"ols"`).

- modality_x, modality_y, channels_x, channels_y:

  Passed to the shared signal-pair extraction when `y` is supplied.

## Value

A list with `dtf` (n x n x n_freqs array), the band-averaged directed
`matrix` (rows = targets, columns = sources), `frequencies`, `channels`,
and `order`. For a two-channel input, `dtf_xy` and `dtf_yx` give the
directional summaries.

## References

Kaminski, M. J., & Blinowska, K. J. (1991). A new method of the
description of the information flow in the brain structures. Biological
Cybernetics, 65(3), 203-210.

## See also

[`partialDirectedCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/partialDirectedCoherence.md),
[`transferEntropy()`](https://x-biosignal.github.io/PhysioCrossModal/reference/transferEntropy.md),
[`grangerCausality()`](https://x-biosignal.github.io/PhysioCrossModal/reference/grangerCausality.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(1500), ncol = 3)
res <- directedTransferFunction(X, sr = 100, order = 4)
dim(res$dtf)
#> [1]   3   3 128
```
