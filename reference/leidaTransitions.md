# LEiDA state-transition metrics

From a
[`leidaStates()`](https://x-biosignal.github.io/PhysioCrossModal/reference/leidaStates.md)
result, computes the state-to-state transition probability matrix, the
switching rate, fractional occupancy, mean dwell time, and the
occupancy-weighted transition entropy (predictability of the next
state).

## Usage

``` r
leidaTransitions(x, sr = NULL)
```

## Arguments

- x:

  A `"leida"` object from
  [`leidaStates()`](https://x-biosignal.github.io/PhysioCrossModal/reference/leidaStates.md).

- sr:

  Optional sampling rate (Hz); if given, dwell times are returned in
  seconds and the switching rate additionally in Hz.

## Value

An object of class `"leida_transitions"`: a list with `transition`
(row-stochastic `state x state` matrix), `counts` (raw transition
counts), `switching_rate` (fraction of steps that change state),
`switching_rate_hz` (switches per second, or `NA`), `occupancy`,
`dwell`, and `entropy` (mean bits of uncertainty about the next state).

## Examples

``` r
set.seed(1); sr <- 100; n <- 2000; t <- seq_len(n) / sr
base <- sin(2 * pi * 10 * t)
X <- cbind(base + rnorm(n, sd = .3), base + rnorm(n, sd = .3),
           sin(2 * pi * 10 * t + 1) + rnorm(n, sd = .3),
           sin(2 * pi * 10 * t + 1) + rnorm(n, sd = .3))
res <- leidaStates(X, sr = sr, freq_band = c(8, 12), n_states = 2, seed = 1)
leidaTransitions(res, sr = sr)
#> <LEiDA state-transition metrics>
#>   states:         2
#>   switching rate: 0.020 per step (2.00 Hz)
#>   transition entropy: 0.14 bits
```
