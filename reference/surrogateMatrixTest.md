# Surrogate-based significance test for coupling matrices

Tests each element of a coupling matrix for significance using surrogate
testing, with correction for multiple comparisons (FDR or Bonferroni).

## Usage

``` r
surrogateMatrixTest(
  x,
  y = NULL,
  method,
  n_surrogates = 199L,
  surrogate_type = c("phase", "timeshift", "iaaft", "aaft"),
  correction = c("fdr", "bonferroni", "none"),
  alpha = 0.05,
  channels_x = NULL,
  channels_y = NULL,
  modality_x = NULL,
  modality_y = NULL,
  cores = 1L,
  ...
)
```

## Arguments

- x:

  PhysioExperiment or MultiPhysioExperiment.

- y:

  PhysioExperiment or NULL (when `x` is an MPE).

- method:

  Character coupling method (same options as
  [`couplingAnalysis`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)).

- n_surrogates:

  Integer number of surrogates (default 199).

- surrogate_type:

  Character surrogate generation method: `"phase"` (default; Fourier
  phase randomization), `"timeshift"` (circular time shift), `"iaaft"`
  (iterative amplitude-adjusted Fourier transform, preserving the power
  spectrum and value distribution; Schreiber & Schmitz, 1996), or
  `"aaft"` (non-iterative amplitude-adjusted FT).

- correction:

  Character correction method: `"fdr"` (default), `"bonferroni"`, or
  `"none"`.

- alpha:

  Numeric significance level (default 0.05).

- channels_x, channels_y:

  Integer vectors of channel indices, or NULL for all channels.

- modality_x, modality_y:

  Character modality names for MPE input.

- cores:

  Integer number of parallel cores (default 1).

- ...:

  Additional arguments passed to the coupling function.

## Value

A list with components:

- matrix:

  Coupling matrix (observed values).

- p_values:

  Matrix of raw p-values (same dimensions).

- p_adjusted:

  Matrix of corrected p-values.

- significant:

  Logical matrix indicating significance.

- correction:

  Character correction method used.

- alpha:

  Numeric significance level.

## References

Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., & Farmer, J. D.
(1992). Testing for nonlinearity in time series: the method of surrogate
data. *Physica D: Nonlinear Phenomena*, 58(1–4), 77–94.

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. *Journal of
the Royal Statistical Society: Series B (Methodological)*, 57(1),
289–300.

## See also

[`surrogateTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateTest.md),
[`couplingMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingMatrix.md),
[`coherenceMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherenceMatrix.md)
