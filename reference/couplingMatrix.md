# Generic coupling matrix across channel pairs

Computes any coupling method for all pairs of channels between two sets
of signals, extracting a scalar statistic for each pair. This is the
multi-channel generalisation of
[`couplingAnalysis`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md).

## Usage

``` r
couplingMatrix(
  x,
  y = NULL,
  sr = NULL,
  method,
  channels_x = NULL,
  channels_y = NULL,
  modality_x = NULL,
  modality_y = NULL,
  ...
)
```

## Arguments

- x:

  A PhysioExperiment or MultiPhysioExperiment object.

- y:

  A PhysioExperiment, or NULL when `x` is an MPE.

- sr:

  Numeric sampling rate in Hz.

- method:

  Character coupling method, one of `"coherence"`, `"plv"`, `"pli"`,
  `"wpli"`, `"granger"`, `"crosscorrelation"`.

- channels_x, channels_y:

  Integer vectors of channel indices, or NULL for all channels.

- modality_x, modality_y:

  Character modality names for MPE input.

- ...:

  Additional arguments passed to the coupling function.

## Value

A list with components:

- matrix:

  Numeric matrix of coupling values.

- method:

  Character method used.

- channel_names_x:

  Character vector of x channel names.

- channel_names_y:

  Character vector of y channel names.

## References

Carter, G. C. (1987). Coherence and time delay estimation. *Proceedings
of the IEEE*, 75(2), 236–255.

## See also

[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md),
[`coherenceMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherenceMatrix.md),
[`plotCouplingMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotCouplingMatrix.md),
[`surrogateMatrixTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateMatrixTest.md)

## Examples

``` r
sr <- 200; n <- sr * 2
pe1 <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
  samplingRate = sr
)
pe2 <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
  samplingRate = sr
)
result <- couplingMatrix(pe1, pe2, method = "coherence", nperseg = 64L)
result$matrix
#>           y_ch1     y_ch2
#> x_ch1 0.2266249 0.3193599
#> x_ch2 0.2407800 0.2505534
```
