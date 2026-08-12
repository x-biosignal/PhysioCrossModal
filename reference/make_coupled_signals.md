# Generate coupled sinusoidal signals

Creates a pair of sinusoidal signals with a shared oscillatory component
at a specified coupling frequency and controllable noise level. Useful
for demonstrating and testing spectral coherence, phase synchrony, and
other cross-modal coupling measures.

## Usage

``` r
make_coupled_signals(
  sr1 = 500,
  sr2 = 500,
  coupling_freq = 20,
  coupling_strength = 0.8,
  noise = 0.2,
  duration = 10
)
```

## Arguments

- sr1:

  Numeric sampling rate in Hz for the first signal (default 500).

- sr2:

  Numeric sampling rate in Hz for the second signal (default 500).

- coupling_freq:

  Numeric frequency in Hz of the shared oscillatory component (default
  20).

- coupling_strength:

  Numeric amplitude of the shared component (default 0.8).

- noise:

  Numeric standard deviation of the additive Gaussian noise (default
  0.2).

- duration:

  Numeric duration in seconds (default 10).

## Value

A named list with components:

- x:

  Numeric vector – first coupled signal.

- y:

  Numeric vector – second coupled signal.

- sr_x:

  Numeric sampling rate of signal x.

- sr_y:

  Numeric sampling rate of signal y.

- coupling_freq:

  Numeric coupling frequency used.

## Details

Both signals share the same sinusoidal component at `coupling_freq`,
scaled by `coupling_strength`, with additive Gaussian noise scaled by
`noise`.

## References

Carter, G. C. (1987). Coherence and time delay estimation. *Proceedings
of the IEEE*, 75(2), 236–255.

## See also

[`make_eeg_emg()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_eeg_emg.md),
[`make_directed_signals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_directed_signals.md),
[`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
signals <- make_coupled_signals(sr1 = 500, sr2 = 500,
                                 coupling_freq = 20,
                                 coupling_strength = 0.8,
                                 noise = 0.2, duration = 10)
result <- coherence(signals$x, signals$y, sr = signals$sr_x)
```
