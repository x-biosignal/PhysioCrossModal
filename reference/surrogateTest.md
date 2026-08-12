# Surrogate-based significance test for coupling

Tests whether the observed coupling between two signals is statistically
significant by comparing it against a null distribution generated from
surrogate signals. Surrogates are created by randomising the phases (or
circularly shifting) one signal, thereby destroying the cross-signal
relationship while preserving the autocorrelation structure.

## Usage

``` r
surrogateTest(
  x,
  y = NULL,
  sr = NULL,
  method,
  n_surrogates = 199L,
  surrogate_type = c("phase", "timeshift", "iaaft", "aaft"),
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L,
  cores = 1L,
  ...
)
```

## Arguments

- x:

  Numeric vector, PhysioExperiment, or MultiPhysioExperiment.

- y:

  Numeric vector or PhysioExperiment, or NULL when `x` is an MPE.

- sr:

  Numeric sampling rate in Hz (required when x/y are numeric).

- method:

  Character coupling method, one of `"coherence"`, `"plv"`, `"pli"`,
  `"wpli"`, `"granger"`, `"crosscorrelation"`, `"wavelet_coherence"`,
  `"wavelet_plv"`.

- n_surrogates:

  Integer number of surrogates (default 199).

- surrogate_type:

  Character surrogate generation method: `"phase"` (default; Fourier
  phase randomization), `"timeshift"` (circular time shift), `"iaaft"`
  (iterative amplitude-adjusted Fourier transform, preserving the power
  spectrum and value distribution; Schreiber & Schmitz, 1996), or
  `"aaft"` (non-iterative amplitude-adjusted FT).

- modality_x, modality_y:

  Character modality names for MPE input.

- channels_x, channels_y:

  Integer channel indices (default 1).

- cores:

  Integer number of parallel cores to use (default `1L`). When
  `cores > 1`, surrogate computations are distributed across cores using
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (Unix) or
  [`parallel::parLapply()`](https://rdrr.io/r/parallel/clusterApply.html)
  (Windows). Falls back to sequential computation on error.

- ...:

  Additional arguments passed to the coupling function (e.g.
  `freq_band`, `nperseg`, `order`).

## Value

A list with components:

- observed:

  The full coupling result from the original signals.

- statistic:

  Numeric scalar: the extracted coupling statistic.

- p_value:

  Numeric scalar: surrogate p-value.

- surrogate_distribution:

  Numeric vector of surrogate statistics.

- threshold_95:

  Numeric scalar: 95th percentile of surrogates.

## Details

The p-value is computed as `(sum(surr >= obs) + 1) / (n_surr + 1)`,
following the conservative correction of Phipson & Smyth (2010).

## References

Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., & Farmer, J. D.
(1992). Testing for nonlinearity in time series: the method of surrogate
data. *Physica D: Nonlinear Phenomena*, 58(1–4), 77–94.

Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be
zero: calculating exact P-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology*, 9(1),
Article 39.

## See also

[`bootstrapCI()`](https://x-biosignal.github.io/PhysioCrossModal/reference/bootstrapCI.md),
[`surrogateMatrixTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateMatrixTest.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
sr <- 500
t <- seq(0, 2, length.out = sr * 2)
x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
y <- 0.8 * sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
result <- surrogateTest(x, y, sr = sr, method = "coherence",
                        n_surrogates = 19, nperseg = 128L)
result$p_value
#> [1] 0.2
```
