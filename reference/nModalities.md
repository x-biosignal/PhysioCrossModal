# Get number of modalities

Get number of modalities

## Usage

``` r
nModalities(x)

# S4 method for class 'MultiPhysioExperiment'
nModalities(x)
```

## Arguments

- x:

  A `MultiPhysioExperiment` object.

## Value

Integer scalar: the number of modalities stored in the object.

## See also

[`modalities()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modalities.md),
[`experiments()`](https://x-biosignal.github.io/PhysioCrossModal/reference/experiments.md)

## Examples

``` r
eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
pe_eeg <- PhysioExperiment(
  assays = list(raw = eeg_data),
  samplingRate = 250
)
mpe <- MultiPhysioExperiment(EEG = pe_eeg)
nModalities(mpe)
#> [1] 1
```
