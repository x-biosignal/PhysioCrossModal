# Show method for MultiPhysioExperiment

Displays a human-readable summary of the object.

## Usage

``` r
# S4 method for class 'MultiPhysioExperiment'
show(object)
```

## Arguments

- object:

  A `MultiPhysioExperiment` object.

## Examples

``` r
eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
pe_eeg <- PhysioExperiment(
  assays = list(raw = eeg_data),
  samplingRate = 250
)
mpe <- MultiPhysioExperiment(EEG = pe_eeg)
mpe
#> class: MultiPhysioExperiment
#> modalities(1): EEG
#> samplingRates: EEG=250Hz
#>   EEG: 500 timepoints x 2 channels
```
