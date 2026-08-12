# Generate simulated EEG-EMG PhysioExperiment pair

Creates a pair of `PhysioExperiment` objects simulating simultaneous EEG
and EMG recordings with corticomuscular coherence (CMC) at a specified
frequency. The first channel of each modality contains a shared
oscillatory component; remaining channels contain independent noise.

## Usage

``` r
make_eeg_emg(
  n_sec = 10,
  n_eeg_ch = 3,
  n_emg_ch = 2,
  sr_eeg = 500,
  sr_emg = 1000,
  cmc_freq = 20
)
```

## Arguments

- n_sec:

  Numeric recording duration in seconds (default 10).

- n_eeg_ch:

  Integer number of EEG channels (default 3).

- n_emg_ch:

  Integer number of EMG channels (default 2).

- sr_eeg:

  Numeric EEG sampling rate in Hz (default 500).

- sr_emg:

  Numeric EMG sampling rate in Hz (default 1000).

- cmc_freq:

  Numeric corticomuscular coherence frequency in Hz (default 20).

## Value

A named list with components:

- eeg:

  A `PhysioExperiment` object with EEG data (`n_eeg_ch` channels at
  `sr_eeg` Hz).

- emg:

  A `PhysioExperiment` object with EMG data (`n_emg_ch` channels at
  `sr_emg` Hz).

## Details

This function is useful for demonstrating and testing multi-channel and
multi-modal coupling analyses such as coherence matrices and
`MultiPhysioExperiment` workflows.

## References

Halliday, D. M., Rosenberg, J. R., Amjad, A. M., Breeze, P., Conway, B.
A., & Farmer, S. F. (1995). A framework for the analysis of mixed time
series/point process data – theory and application to the study of
physiological tremor, single motor unit discharges and electromyograms.
*Progress in Biophysics and Molecular Biology*, 64(2–3), 237–278.

## See also

[`make_coupled_signals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_coupled_signals.md),
[`make_directed_signals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_directed_signals.md),
[`MultiPhysioExperiment()`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md),
[`coherenceMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherenceMatrix.md)

## Examples

``` r
data <- make_eeg_emg(n_sec = 5, n_eeg_ch = 3, n_emg_ch = 2)
mpe <- MultiPhysioExperiment(EEG = data$eeg, EMG = data$emg)
modalities(mpe)
#> [1] "EEG" "EMG"
```
