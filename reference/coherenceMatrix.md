# Coherence matrix across channel pairs

Computes magnitude-squared coherence for all pairs of channels between
two sets of signals (or between all channels within a single set).
Returns a matrix of peak coherence values plus the full spectra.

## Usage

``` r
coherenceMatrix(
  x,
  y = NULL,
  sr = NULL,
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

  Numeric sampling rate in Hz (used only for numeric inputs).

- channels_x:

  Integer vector of channel indices for x, or NULL for all channels.

- channels_y:

  Integer vector of channel indices for y, or NULL for all channels.

- modality_x, modality_y:

  Character modality names for MPE input.

- ...:

  Additional arguments passed to
  [`coherence`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md)
  (e.g. `nperseg`, `freq_range`).

## Value

A list with components:

- matrix:

  Numeric matrix of peak coherence values, with dimensions
  `[n_channels_x x n_channels_y]`.

- spectra:

  List of lists containing the full coherence result for each pair.

- frequencies:

  Numeric vector of frequencies from the first pair.

- channel_names_x:

  Character vector of x channel names.

- channel_names_y:

  Character vector of y channel names.

## References

Carter, G. C. (1987). Coherence and time delay estimation. *Proceedings
of the IEEE*, 75(2), 236–255.

## See also

[`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md),
[`couplingMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingMatrix.md),
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
result <- coherenceMatrix(pe1, pe2, nperseg = 64L)
result$matrix
#>           y_ch1     y_ch2
#> x_ch1 0.3469981 0.4737342
#> x_ch2 0.3909215 0.3081293
```
