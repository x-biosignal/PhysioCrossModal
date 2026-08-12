# Granger Causality Between Two Signals

Computes directed coupling between two physiological signals using
either autoregressive prediction or nonparametric spectral
factorization. Both scalar and frequency-resolved Granger causality are
supported.

## Usage

``` r
grangerCausality(
  x,
  y = NULL,
  sr = NULL,
  order = 5L,
  freq_range = NULL,
  method = c("parametric", "nonparametric"),
  modality_x = NULL,
  modality_y = NULL,
  channels_x = 1L,
  channels_y = 1L,
  spectral_estimator = c("multitaper", "welch"),
  n_fft = NULL,
  segment_length = NULL,
  overlap = 0.5,
  time_bandwidth = 4,
  n_tapers = NULL,
  factor_tolerance = 1e-08,
  factor_max_iterations = 1000L,
  eigen_floor = 1e-10
)
```

## Arguments

- x:

  Numeric vector, PhysioExperiment, or MultiPhysioExperiment.

- y:

  Numeric vector or PhysioExperiment, or NULL when `x` is a
  MultiPhysioExperiment.

- sr:

  Numeric sampling rate in Hz (required when x/y are numeric vectors).

- order:

  Integer AR model order (default 5). Higher order captures longer-range
  temporal dependencies but requires more data.

- freq_range:

  Numeric vector of length 2 specifying the frequency band in Hz over
  which to compute spectral Granger causality. With the parametric
  method, `NULL` returns time-domain GC. With the nonparametric method,
  `NULL` aggregates the full DC-to-Nyquist spectrum.

- method:

  Character estimator. `"parametric"` (default) uses the established
  MVAR ordinary-least-squares path. `"nonparametric"` estimates a
  complete cross-spectral matrix and applies Wilson factorization
  followed by Geweke's frequency-domain decomposition.

- modality_x:

  Character modality name in MultiPhysioExperiment for the x signal.

- modality_y:

  Character modality name in MultiPhysioExperiment for the y signal.

- channels_x:

  Integer channel index to extract from x (default 1).

- channels_y:

  Integer channel index to extract from y (default 1).

- spectral_estimator:

  Nonparametric cross-spectral estimator: `"multitaper"` (default) or
  `"welch"`.

- n_fft:

  Even FFT length for the nonparametric estimator. `NULL` chooses a
  deterministic length from `segment_length`.

- segment_length:

  Segment length for cross-spectral estimation. `NULL` chooses the
  largest usable power of two up to 512 samples.

- overlap:

  Fractional segment overlap in `[0, 1)`.

- time_bandwidth:

  DPSS time-bandwidth product for multitaper estimation.

- n_tapers:

  Number of DPSS tapers; `NULL` uses `floor(2 * time_bandwidth) - 1`.

- factor_tolerance:

  Relative Wilson-factor update tolerance.

- factor_max_iterations:

  Maximum Wilson iterations.

- eigen_floor:

  Relative positive eigenvalue floor used only for round-off-scale or
  ill-conditioned cross-spectral bins.

## Value

A list with components:

- gc_xy:

  Numeric Granger causality from x to y (x drives y). A positive value
  indicates that x's past helps predict y beyond y's own past.

- gc_yx:

  Numeric Granger causality from y to x (y drives x).

- net_gc:

  Numeric net Granger causality (`gc_xy - gc_yx`). Positive values
  indicate that x drives y more than y drives x.

- order:

  Integer AR model order used, or `NA` for nonparametric estimation.

- method:

  Character estimator used.

Nonparametric results additionally contain the frequency grid,
directional spectra, estimator/factorization metadata, innovation
covariance, the factorization reconstruction limit, and regularization
diagnostics.

## Details

For the time-domain (parametric) method, Granger causality from x to y
is defined as: \$\$GC\_{x \to y} = \log(\sigma^2\_{restricted} /
\sigma^2\_{unrestricted})\$\$ where the restricted model predicts y from
its own past only, and the unrestricted model predicts y from both its
own past and x's past.

For spectral Granger causality (when `freq_range` is specified), the AR
coefficients are transformed to the frequency domain via the transfer
function \\H(f) = A(f)^{-1}\\, and the spectral GC is computed as:
\$\$GC\_{x \to y}(f) = \log(S\_{yy}^{restricted}(f) /
S\_{yy}^{unrestricted}(f))\$\$ The returned values are then averaged
over the specified frequency band.

The parametric frequency-band path uses separate restricted and
unrestricted univariate AR models and is an approximation to the full
Geweke spectral decomposition. Its established output is retained for
compatibility. The nonparametric path instead estimates a complete
Hermitian 2x2 cross-spectral matrix with common windows/tapers, applies
Wilson's minimum-phase factorization, and evaluates Geweke directional
spectra. A scalar is the trapezoidal frequency mean over the selected
band, with DC and Nyquist receiving endpoint half weights.

Cross-spectral eigenvalue adjustment is limited to the reported relative
`eigen_floor`; materially indefinite matrices fail. Wilson update
convergence and spectral reconstruction are checked separately. The
relative reconstruction limit is
`min(0.05, max(0.02, 100 * factor_tolerance))`, so it remains between 2%
and 5%; the observed error and limit are returned in the result.
Nonparametric estimation requires at least 32 paired samples and at
least two spectral realizations.

Granger causality is a stationary linear predictive-direction measure.
It is not proof of physical, mechanistic, or interventional causality
and does not resolve unmeasured confounding, common drivers, volume
conduction, instantaneous mixing, or non-stationarity.

## References

Granger, C. W. J. (1969). Investigating causal relations by econometric
models and cross-spectral methods. Econometrica, 37(3), 424-438.

Geweke, J. (1982). Measurement of linear dependence and feedback between
multiple time series. Journal of the American Statistical Association,
77(378), 304-313.

Wilson, G. T. (1972). The factorization of matricial spectral densities.
SIAM Journal on Applied Mathematics, 23(4), 420-426.

Dhamala, M., Rangarajan, G., & Ding, M. (2008). Analyzing information
flow in brain networks with nonparametric Granger causality. NeuroImage,
41(2), 354-362.

Dhamala, M., Rangarajan, G., & Ding, M. (2008). Estimating Granger
causality from Fourier and wavelet transforms of time series data.
Physical Review Letters, 100, 018701.

## See also

[`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md),
[`phaseLockingValue()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLockingValue.md)

## Examples

``` r
# Create directed signals: x drives y with a 5-sample lag
set.seed(42)
n <- 5000
x <- rnorm(n)
y <- numeric(n)
for (i in 6:n) y[i] <- 0.6 * x[i - 5] + 0.4 * rnorm(1)

result <- grangerCausality(x, y, sr = 500, order = 10)
result$gc_xy   # Should be positive (x drives y)
#> [1] 1.173665
result$net_gc  # Should be positive
#> [1] 1.170759
```
