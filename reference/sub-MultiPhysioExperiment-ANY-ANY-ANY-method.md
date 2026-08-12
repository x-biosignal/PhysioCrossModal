# Subset a MultiPhysioExperiment by time range and/or modality

Provides flexible subsetting of `MultiPhysioExperiment` objects:

- `mpe[, "eeg"]` – select modalities by name, returns a new MPE

- `mpe[c(1.0, 5.0), ]` – select a time window in seconds across all
  modalities, accounting for each modality's sampling rate

- `mpe[c(1.0, 5.0), "eeg"]` – both time and modality subsetting

## Usage

``` r
# S4 method for class 'MultiPhysioExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  A `MultiPhysioExperiment` object.

- i:

  Numeric vector of length 2 specifying time range `[tmin, tmax]` in
  seconds, or missing.

- j:

  Character vector of modality names, or missing.

- ...:

  Additional arguments (not used).

- drop:

  Logical (not used).

## Value

A `MultiPhysioExperiment` object with the selected modalities and/or
time windows.

## Details

When `i` is a numeric vector of length 2, it is interpreted as a time
range `[tmin, tmax]` in seconds. For each modality, the appropriate
sample indices are computed from its sampling rate:
`start_idx = floor(tmin * sr) + 1`, `end_idx = floor(tmax * sr) + 1`,
clamped to valid bounds.

## Examples

``` r
eeg_data <- matrix(rnorm(2500 * 4), nrow = 2500, ncol = 4)
emg_data <- matrix(rnorm(5000 * 2), nrow = 5000, ncol = 2)
pe_eeg <- PhysioExperiment(
  assays = list(raw = eeg_data),
  samplingRate = 500
)
pe_emg <- PhysioExperiment(
  assays = list(raw = emg_data),
  samplingRate = 1000
)
mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)

# Subset by modality
mpe_eeg <- mpe[, "EEG"]
modalities(mpe_eeg)  # "EEG"
#> [1] "EEG"

# Subset by time range (seconds)
mpe_win <- mpe[c(1.0, 3.0), ]

# Combined
mpe_sub <- mpe[c(1.0, 3.0), "EEG"]
```
