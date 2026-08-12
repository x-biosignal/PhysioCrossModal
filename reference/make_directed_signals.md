# Generate directed coupling signals

Creates a pair of signals where `x` drives `y` with a specified lag and
coupling strength. Signal `x` is white noise, and `y` is a mixture of a
lagged copy of `x` and independent noise. This is useful for testing
directed coupling measures such as Granger causality.

## Usage

``` r
make_directed_signals(n = 5000, sr = 500, lag_samples = 10, coupling = 0.7)
```

## Arguments

- n:

  Integer number of samples (default 5000).

- sr:

  Numeric sampling rate in Hz (default 500).

- lag_samples:

  Integer number of samples by which `x` leads `y` (default 10).

- coupling:

  Numeric coupling strength in \[0, 1\] (default 0.7).

## Value

A named list with components:

- x:

  Numeric vector – driving signal (white noise).

- y:

  Numeric vector – driven signal (lagged mixture).

- sr:

  Numeric sampling rate.

- lag_samples:

  Integer lag used.

## Details

The generating model is: \$\$y(t) = \text{coupling} \cdot x(t -
\text{lag\\samples}) + (1 - \text{coupling}) \cdot \epsilon(t)\$\$ where
\\\epsilon(t) \sim N(0, 1)\\.

## References

Granger, C. W. J. (1969). Investigating causal relations by econometric
models and cross-spectral methods. *Econometrica*, 37(3), 424–438.

## See also

[`grangerCausality()`](https://x-biosignal.github.io/PhysioCrossModal/reference/grangerCausality.md),
[`make_coupled_signals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_coupled_signals.md),
[`crossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossCorrelation.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
signals <- make_directed_signals(n = 5000, sr = 500,
                                  lag_samples = 10, coupling = 0.7)
result <- grangerCausality(signals$x, signals$y, sr = signals$sr,
                           order = 15)
result$gc_xy   # should be positive (x drives y)
#> [1] 1.846288
result$net_gc  # should be positive
#> [1] 1.843261
```
