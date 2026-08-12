# Get sampling rates for all modalities

Returns a named numeric vector of sampling rates, one per modality.

## Usage

``` r
samplingRates(x)

# S4 method for class 'MultiPhysioExperiment'
samplingRates(x)
```

## Arguments

- x:

  A `MultiPhysioExperiment` object.

## Value

Named numeric vector of sampling rates in Hz (e.g.,
`c(EEG = 500, EMG = 1000)`).

## See also

[`modalities()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modalities.md),
[`experiments()`](https://x-biosignal.github.io/PhysioCrossModal/reference/experiments.md),
[`alignToRate()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignToRate.md)

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
samplingRates(mpe)
#>  EEG  EMG 
#>  250 1000 
```
