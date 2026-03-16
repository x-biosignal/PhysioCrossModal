# Synthetic data generators for PhysioCrossModal tests

#' Create two coupled sinusoidal signals at different sampling rates
make_coupled_signals <- function(sr1 = 500, sr2 = 1000, n_sec = 10,
                                  coupling_freq = 20, coupling_strength = 0.8,
                                  noise_level = 0.2, seed = 42) {
  set.seed(seed)
  n1 <- as.integer(sr1 * n_sec)
  n2 <- as.integer(sr2 * n_sec)
  t1 <- seq(0, n_sec, length.out = n1)
  t2 <- seq(0, n_sec, length.out = n2)

  driver1 <- sin(2 * pi * coupling_freq * t1)
  driver2 <- sin(2 * pi * coupling_freq * t2)

  x <- coupling_strength * driver1 + noise_level * rnorm(n1)
  y <- coupling_strength * driver2 + noise_level * rnorm(n2)

  list(x = x, y = y, sr_x = sr1, sr_y = sr2, t_x = t1, t_y = t2)
}

#' Create PhysioExperiment pair simulating EEG and EMG
make_eeg_emg <- function(n_sec = 10, n_eeg_ch = 4, n_emg_ch = 2,
                          sr_eeg = 500, sr_emg = 1000,
                          cmc_freq = 20, cmc_strength = 0.7,
                          seed = 42) {
  set.seed(seed)
  n_eeg <- as.integer(sr_eeg * n_sec)
  n_emg <- as.integer(sr_emg * n_sec)
  t_eeg <- seq(0, n_sec, length.out = n_eeg)
  t_emg <- seq(0, n_sec, length.out = n_emg)

  driver_eeg <- sin(2 * pi * cmc_freq * t_eeg)
  driver_emg <- sin(2 * pi * cmc_freq * t_emg)

  eeg_data <- matrix(rnorm(n_eeg * n_eeg_ch), nrow = n_eeg, ncol = n_eeg_ch)
  eeg_data[, 1] <- eeg_data[, 1] + cmc_strength * driver_eeg

  emg_data <- matrix(rnorm(n_emg * n_emg_ch) * 0.5, nrow = n_emg, ncol = n_emg_ch)
  emg_data[, 1] <- emg_data[, 1] + cmc_strength * driver_emg

  pe_eeg <- PhysioExperiment(
    assays = list(raw = eeg_data),
    colData = S4Vectors::DataFrame(
      label = paste0("EEG", seq_len(n_eeg_ch)),
      type = rep("EEG", n_eeg_ch)
    ),
    samplingRate = sr_eeg
  )

  pe_emg <- PhysioExperiment(
    assays = list(raw = emg_data),
    colData = S4Vectors::DataFrame(
      label = paste0("EMG", seq_len(n_emg_ch)),
      type = rep("EMG", n_emg_ch)
    ),
    samplingRate = sr_emg
  )

  list(eeg = pe_eeg, emg = pe_emg)
}

#' Create signals with known directional coupling (x drives y with lag)
make_directed_signals <- function(n = 5000, sr = 500, lag_samples = 5,
                                   coupling = 0.6, seed = 42) {
  set.seed(seed)
  x <- rnorm(n)
  y <- numeric(n)
  for (i in (lag_samples + 1):n) {
    y[i] <- coupling * x[i - lag_samples] + (1 - coupling) * rnorm(1)
  }
  list(x = x, y = y, sr = sr, lag_samples = lag_samples)
}
