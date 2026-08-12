# Unified interface for cross-modal coupling analysis

Dispatches to the appropriate coupling function based on the `method`
parameter. Accepts the same flexible inputs as the underlying functions:
two numeric vectors (with `sr`), two PhysioExperiment objects, or a
MultiPhysioExperiment with named modalities.

## Usage

``` r
couplingAnalysis(
  x,
  y = NULL,
  mpe = NULL,
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L,
  method = c("coherence", "plv", "pli", "wpli", "granger", "crosscorrelation",
    "wavelet_coherence", "wavelet_plv", "multitaper_coherence", "dtf", "pdc",
    "transferentropy", "pac"),
  sr = NULL,
  ...,
  granger_method = c("parametric", "nonparametric")
)
```

## Arguments

- x:

  Numeric vector, PhysioExperiment, or MultiPhysioExperiment. When `mpe`
  is provided this argument is ignored.

- y:

  Numeric vector or PhysioExperiment (NULL when `x` is an MPE or when
  `mpe` is provided).

- mpe:

  A
  [`MultiPhysioExperiment`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md)
  object. When supplied, `x` and `y` are ignored and signals are
  extracted from `mpe` using `modality_x` / `modality_y`.

- modality_x, modality_y:

  Character names of the modalities to extract from `mpe`.

- channels_x, channels_y:

  Integer channel indices to extract (default 1).

- method:

  Character string specifying the coupling method. One of:
  `"coherence"`, `"plv"`, `"pli"`, `"wpli"`, `"granger"`,
  `"crosscorrelation"`, `"wavelet_coherence"`, or `"wavelet_plv"`.

- sr:

  Numeric sampling rate in Hz. Required when `x` and `y` are numeric
  vectors.

- ...:

  Additional arguments passed to the specific coupling function (e.g.
  `freq_band`, `order`, `max_lag`, `nperseg`, etc.).

- granger_method:

  Named-only estimator selector used when `method = "granger"`:
  `"parametric"` (default) or `"nonparametric"`.

## Value

The result from the dispatched coupling function. See individual
function documentation for details:

- `"coherence"`: see
  [`coherence`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md)

- `"plv"`: see
  [`phaseLockingValue`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLockingValue.md)

- `"pli"`: see
  [`phaseLagIndex`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLagIndex.md)

- `"wpli"`: see
  [`weightedPLI`](https://x-biosignal.github.io/PhysioCrossModal/reference/weightedPLI.md)

- `"granger"`: see
  [`grangerCausality`](https://x-biosignal.github.io/PhysioCrossModal/reference/grangerCausality.md)

- `"crosscorrelation"`: see
  [`crossCorrelation`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossCorrelation.md)

- `"wavelet_coherence"`: see
  [`waveletCoherence`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletCoherence.md)

- `"wavelet_plv"`: see
  [`waveletPLV`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletPLV.md)

## References

Carter, G. C. (1987). Coherence and time delay estimation. *Proceedings
of the IEEE*, 75(2), 236–255.

Lachaux, J.-P., Rodriguez, E., Martinerie, J., & Varela, F. J. (1999).
Measuring phase synchrony in brain signals. *Human Brain Mapping*, 8(4),
194–208.

Granger, C. W. J. (1969). Investigating causal relations by econometric
models and cross-spectral methods. *Econometrica*, 37(3), 424–438.

## See also

[`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md),
[`phaseLockingValue()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLockingValue.md),
[`grangerCausality()`](https://x-biosignal.github.io/PhysioCrossModal/reference/grangerCausality.md),
[`crossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossCorrelation.md),
[`surrogateTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateTest.md)

## Examples

``` r
# Numeric vectors
sr <- 500
t <- seq(0, 10, length.out = sr * 10)
x <- sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
y <- 0.8 * sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))

# Coherence
res <- couplingAnalysis(x, y, method = "coherence", sr = sr)

# Cross-correlation
res <- couplingAnalysis(x, y, method = "crosscorrelation", sr = sr)
```
