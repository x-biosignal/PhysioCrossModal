# Align multiple PhysioExperiment objects to a common sampling rate

Takes two or more named `PhysioExperiment` objects and resamples them so
they all share a single sampling rate, then wraps the result in a
[`MultiPhysioExperiment`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md).

## Usage

``` r
alignSignals(
  ...,
  method = c("lowest_rate", "common_rate", "resample", "elastic"),
  target_rate = NULL
)
```

## Arguments

- ...:

  Named `PhysioExperiment` objects.

- method:

  Character: one of `"lowest_rate"` (default), `"common_rate"`,
  `"resample"`, or `"elastic"`.

- target_rate:

  Numeric scalar: required when `method = "resample"`.

## Value

A
[`MultiPhysioExperiment`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md)
containing the (possibly resampled) input objects.

## Details

Three alignment strategies are available:

- `"lowest_rate"`:

  (default) Resample all signals to the lowest sampling rate present.

- `"common_rate"`:

  Resample all signals to the highest sampling rate present.

- `"resample"`:

  Resample all signals to the rate given by `target_rate` (which must be
  supplied).

- `"elastic"`:

  Resample to the lowest rate, then additionally apply an elastic (SRVF)
  time warping of every object onto the first, so that shared events
  line up in phase. The warping is computed from each object's first
  channel and carried to all its channels; it is stored in
  `metadata()$elastic_warp`. Objects whose duration differs from the
  reference after rate alignment are left un-warped. See
  [`elasticAlign()`](https://x-biosignal.github.io/PhysioCrossModal/reference/elasticAlign.md).

If every input already shares the same sampling rate, no resampling is
performed regardless of the chosen method (elastic warping, if
requested, is still applied).

## See also

[`alignToRate()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignToRate.md),
[`mergePhysio()`](https://x-biosignal.github.io/PhysioCrossModal/reference/mergePhysio.md),
[`MultiPhysioExperiment()`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
eeg <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(500 * 2), nrow = 500)),
  samplingRate = 500
)
emg <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000 * 2), nrow = 1000)),
  samplingRate = 1000
)
mpe <- alignSignals(EEG = eeg, EMG = emg, method = "lowest_rate")
```
