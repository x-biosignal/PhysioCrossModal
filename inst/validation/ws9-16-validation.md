# WS9-16 independent numeric validation

Pinned external fixture: Tensorpac 0.6.5 (Tort and Canolty); Ozkurt uses
the independently coded published normalization because Tensorpac's
`norm_direct_pac` is a different statistic.

| gate | value | threshold | status |
|---|---:|---:|---|
| fixture_hashes | 1e+00 | == 1 | PASS |
| external_formula_max_abs_error | 3.330669074e-16 | <= 1e-12 | PASS |
| independent_formula_max_abs_error | 2.775557562e-17 | <= 1e-12 | PASS |
| public_grid_max_abs_error | 5.551115123e-17 | <= 1e-10 | PASS |
| coupled_adjusted_sensitivity | 1e+00 | >= 0.80 | PASS |
| uncoupled_adjusted_cell_fpr | 5e-03 | <= 0.10 | PASS |
| scale_invariance_max_error | 8.881784197e-16 | <= 1e-12 | PASS |
| rng_restored | 1e+00 | == 1 | PASS |
| mutation_detection | 7e+00 | == 7 | PASS |

- Formula fixtures: 100
- Coupled simulations: 100
- Uncoupled simulations: 100 (400 tested cells)
- Surrogates per simulation: 199
- Mutations detected: 7/7
- Quantile/p-value convention: upper tail with Phipson-Smyth +1 correction.
- Filter reference: signal 1.8.1 Butterworth + filtfilt; Hilbert transform
  independently reconstructed from the FFT analytic-signal multiplier.
