# LEiDA dynamic functional-connectivity states

Leading Eigenvector Dynamics Analysis (Cabral et al. 2017): at each time
point the instantaneous phase-coherence matrix is reduced to its leading
eigenvector; these are clustered over time into a small set of recurrent
connectivity states. Complements state-resolved network analysis by
finding the states data-drivenly rather than from external labels.

## Usage

``` r
leidaStates(
  x,
  sr,
  freq_band,
  n_states = 4L,
  order = 4L,
  nstart = 10L,
  seed = NULL
)
```

## Arguments

- x:

  A numeric matrix, time (rows) by channels (columns).

- sr:

  Sampling rate in Hz.

- freq_band:

  Numeric `c(low, high)` band in Hz.

- n_states:

  Number of clusters/states (default 4).

- order:

  Butterworth filter order (default 4).

- nstart:

  k-means restarts (default 10).

- seed:

  Optional integer seed for reproducible clustering.

## Value

An object of class `"leida"`: a list with `states` (per-time cluster
labels), `centroids` (states x channels leading-eigenvector centroids),
`occupancy` (fraction of time in each state), `dwell` (mean run length
in samples per state), and `n_states`.

## References

Cabral J, et al. (2017). Cognitive performance in healthy older adults
relates to spontaneous switching between states of functional
connectivity during rest. *Sci Rep* 7:5135.

## Examples

``` r
set.seed(1); sr <- 100; n <- 2000; t <- seq_len(n) / sr
# two regimes: channels 1-2 coherent first half, 3-4 coherent second half
base <- sin(2 * pi * 10 * t)
X <- cbind(base + rnorm(n, sd = .3), base + rnorm(n, sd = .3),
           sin(2 * pi * 10 * t + 1) + rnorm(n, sd = .3),
           sin(2 * pi * 10 * t + 1) + rnorm(n, sd = .3))
res <- leidaStates(X, sr = sr, freq_band = c(8, 12), n_states = 2, seed = 1)
res$occupancy
#> [1] 0.43 0.57
```
