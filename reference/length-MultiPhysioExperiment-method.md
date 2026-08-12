# Additional methods for MultiPhysioExperiment

Provides subsetting (`[`),
[`length()`](https://rdrr.io/r/base/length.html), and
[`names()`](https://rdrr.io/r/base/names.html) methods for
[`MultiPhysioExperiment`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md)
objects. Length of a MultiPhysioExperiment

## Usage

``` r
# S4 method for class 'MultiPhysioExperiment'
length(x)
```

## Arguments

- x:

  A `MultiPhysioExperiment` object.

## Value

Integer: number of modalities.

## Details

Returns the number of modalities (same as
[`nModalities`](https://x-biosignal.github.io/PhysioCrossModal/reference/nModalities.md)).

## Examples

``` r
eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
pe_eeg <- PhysioExperiment(
  assays = list(raw = eeg_data),
  samplingRate = 250
)
mpe <- MultiPhysioExperiment(EEG = pe_eeg)
length(mpe)  # 1
#> [1] 1
```
