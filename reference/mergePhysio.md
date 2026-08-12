# Merge two PhysioExperiment objects by combining channels

Horizontally concatenates the channels (columns) of two
`PhysioExperiment` objects that share the **same** sampling rate.
Channel labels are prefixed to avoid name collisions.

## Usage

``` r
mergePhysio(x, y, prefix = c("x_", "y_"))
```

## Arguments

- x:

  A `PhysioExperiment` object.

- y:

  A `PhysioExperiment` object.

- prefix:

  Character vector of length 2: prefixes added to channel labels from
  `x` and `y` respectively (default `c("x_", "y_")`).

## Value

A single `PhysioExperiment` with the columns of both inputs.

## Details

Both objects must have the same number of time points (rows) and the
same sampling rate. Only the first assay of each object is merged.

## See also

[`alignToRate()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignToRate.md),
[`alignSignals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignSignals.md),
[`MultiPhysioExperiment()`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md)

## Examples

``` r
pe1 <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(100 * 4), nrow = 100)),
  colData = S4Vectors::DataFrame(label = paste0("A", 1:4)),
  samplingRate = 500
)
pe2 <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(100 * 4), nrow = 100)),
  colData = S4Vectors::DataFrame(label = paste0("B", 1:4)),
  samplingRate = 500
)
merged <- mergePhysio(pe1, pe2)
```
