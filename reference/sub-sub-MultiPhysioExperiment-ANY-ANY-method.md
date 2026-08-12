# Extract a single modality from a MultiPhysioExperiment

Extract a single modality from a MultiPhysioExperiment

## Usage

``` r
# S4 method for class 'MultiPhysioExperiment,ANY,ANY'
x[[i, j, ...]]
```

## Arguments

- x:

  A `MultiPhysioExperiment` object.

- i:

  Modality name (character) or index (integer).

- j:

  Not used.

- ...:

  Additional arguments (not used).

## Value

A `PhysioExperiment` object.

## Examples

``` r
eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
pe_eeg <- PhysioExperiment(
  assays = list(raw = eeg_data),
  samplingRate = 250
)
mpe <- MultiPhysioExperiment(EEG = pe_eeg)
mpe[["EEG"]]
#> class: PhysioExperiment
#> dim: 500 x 2 
#> assays(1): raw
#> samplingRate: 250 Hz
#> channels(2): Ch1, Ch2
```
