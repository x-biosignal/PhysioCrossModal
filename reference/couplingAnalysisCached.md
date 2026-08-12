# Cached coupling analysis over a MultiPhysioExperiment

Wraps
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)
with a per-object memoisation cache held in the `couplingResults` slot
of a
[`MultiPhysioExperiment`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md).
Results are keyed by a
[`digest`](https://eddelbuettel.github.io/digest/man/digest.html) of the
analysis arguments (method, modalities, channels, and any extra
arguments), so a repeated identical request returns the stored result
without recomputing.

## Usage

``` r
couplingAnalysisCached(
  mpe,
  method,
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L,
  ...,
  cache = TRUE
)
```

## Arguments

- mpe:

  A `MultiPhysioExperiment` object.

- method:

  Character coupling method (see
  [`couplingAnalysis`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)).

- modality_x, modality_y:

  Character modality names to extract from `mpe`.

- channels_x, channels_y:

  Integer channel indices (default 1).

- ...:

  Additional arguments forwarded to
  [`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)
  and folded into the cache key.

- cache:

  Logical; if `TRUE` (default) reuse and populate the cache. `FALSE`
  forces recomputation but still stores the fresh result.

## Value

A list with `result` (the coupling result), `mpe` (the object with its
`couplingResults` cache updated) and `cached` (`TRUE` when the result
was served from the cache).

## See also

[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md),
[`couplingResults()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingResults.md)

## Examples

``` r
eeg <- PhysioExperiment(assays = list(raw = matrix(rnorm(1000), 500, 2)),
                        samplingRate = 250)
emg <- PhysioExperiment(assays = list(raw = matrix(rnorm(1000), 500, 2)),
                        samplingRate = 250)
mpe <- MultiPhysioExperiment(EEG = eeg, EMG = emg)
out  <- couplingAnalysisCached(mpe, "coherence", "EEG", "EMG")
out2 <- couplingAnalysisCached(out$mpe, "coherence", "EEG", "EMG")
out2$cached  # TRUE
#> [1] TRUE
```
