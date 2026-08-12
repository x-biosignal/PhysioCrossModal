# Get or set temporal alignment metadata

Get or set temporal alignment metadata

## Usage

``` r
alignment(x)

# S4 method for class 'MultiPhysioExperiment'
alignment(x)

alignment(x) <- value

# S4 method for class 'MultiPhysioExperiment'
alignment(x) <- value
```

## Arguments

- x:

  A `MultiPhysioExperiment` object.

- value:

  A
  [`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  with alignment metadata.

## Value

A [`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html).

## See also

[`alignSignals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignSignals.md),
[`experiments()`](https://x-biosignal.github.io/PhysioCrossModal/reference/experiments.md),
[`samplingRates()`](https://x-biosignal.github.io/PhysioCrossModal/reference/samplingRates.md)

## Examples

``` r
eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
pe_eeg <- PhysioExperiment(
  assays = list(raw = eeg_data),
  samplingRate = 250
)
mpe <- MultiPhysioExperiment(EEG = pe_eeg)
alignment(mpe)
#> DataFrame with 1 row and 3 columns
#>      modality samplingRate    offset
#>   <character>    <numeric> <numeric>
#> 1         EEG          250         0
```
