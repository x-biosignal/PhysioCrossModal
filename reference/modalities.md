# Get modality names

Returns a character vector of the names of the modalities stored in the
object.

## Usage

``` r
modalities(x)

# S4 method for class 'MultiPhysioExperiment'
modalities(x)
```

## Arguments

- x:

  A `MultiPhysioExperiment` object.

## Value

Character vector of modality names (e.g., `c("EEG", "EMG")`).

## See also

[`experiments()`](https://x-biosignal.github.io/PhysioCrossModal/reference/experiments.md),
[`nModalities()`](https://x-biosignal.github.io/PhysioCrossModal/reference/nModalities.md),
[`samplingRates()`](https://x-biosignal.github.io/PhysioCrossModal/reference/samplingRates.md)

## Examples

``` r
eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
emg_data <- matrix(rnorm(1000 * 2), nrow = 1000, ncol = 2)
pe_eeg <- PhysioExperiment(
  assays = list(raw = eeg_data),
  samplingRate = 250
)
pe_emg <- PhysioExperiment(
  assays = list(raw = emg_data),
  samplingRate = 1000
)
mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)
modalities(mpe)
#> [1] "EEG" "EMG"
```
