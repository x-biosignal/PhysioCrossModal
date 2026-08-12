# Get or set the experiments list

Get or set the experiments list

## Usage

``` r
experiments(x)

# S4 method for class 'MultiPhysioExperiment'
experiments(x)

experiments(x) <- value

# S4 method for class 'MultiPhysioExperiment'
experiments(x) <- value
```

## Arguments

- x:

  A `MultiPhysioExperiment` object.

- value:

  A named list of `PhysioExperiment` objects.

## Value

For the getter, a named list of `PhysioExperiment` objects. For the
setter, the modified `MultiPhysioExperiment` (returned invisibly).

## See also

[`MultiPhysioExperiment()`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md),
[`modalities()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modalities.md),
[`samplingRates()`](https://x-biosignal.github.io/PhysioCrossModal/reference/samplingRates.md)

## Examples

``` r
eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
pe_eeg <- PhysioExperiment(
  assays = list(raw = eeg_data),
  samplingRate = 250
)
mpe <- MultiPhysioExperiment(EEG = pe_eeg)
experiments(mpe)
#> $EEG
#> class: PhysioExperiment
#> dim: 500 x 2 
#> assays(1): raw
#> samplingRate: 250 Hz
#> channels(2): Ch1, Ch2
#> 
```
