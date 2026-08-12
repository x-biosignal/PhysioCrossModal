# Resample a PhysioExperiment to a target sampling rate

Resamples all assays in a `PhysioExperiment` object so that the
resulting data correspond to the specified target sampling rate. If the
PhysioPreprocess package is available and `method = "linear"`, its
`resample()` function is used; otherwise an internal interpolation is
applied.

## Usage

``` r
alignToRate(x, target_rate, method = c("linear", "spline", "fft"))
```

## Arguments

- x:

  A `PhysioExperiment` object.

- target_rate:

  Numeric scalar: desired sampling rate in Hz.

- method:

  Character: resampling method. One of `"linear"` (default), `"spline"`,
  or `"fft"`.

## Value

A new `PhysioExperiment` with resampled data and updated `samplingRate`.
The assay name is preserved from the original.

## See also

[`alignSignals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignSignals.md),
[`mergePhysio()`](https://x-biosignal.github.io/PhysioCrossModal/reference/mergePhysio.md),
[`MultiPhysioExperiment()`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md)

## Examples

``` r
eeg_data <- matrix(rnorm(1000 * 4), nrow = 1000, ncol = 4)
pe <- PhysioExperiment(
  assays = list(raw = eeg_data),
  colData = S4Vectors::DataFrame(
    label = paste0("Ch", 1:4),
    type = rep("EEG", 4)
  ),
  samplingRate = 1000
)
pe500 <- alignToRate(pe, target_rate = 500)
samplingRate(pe500)  # 500
#> [1] 500
```
