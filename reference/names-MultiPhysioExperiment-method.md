# Names of a MultiPhysioExperiment

Returns modality names (same as
[`modalities`](https://x-biosignal.github.io/PhysioCrossModal/reference/modalities.md)).

## Usage

``` r
# S4 method for class 'MultiPhysioExperiment'
names(x)
```

## Arguments

- x:

  A `MultiPhysioExperiment` object.

## Value

Character vector of modality names.

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
names(mpe)  # c("EEG", "EMG")
#> [1] "EEG" "EMG"
```
