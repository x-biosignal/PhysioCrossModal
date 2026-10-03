# Exported synthetic data generators for PhysioCrossModal
#
# Convenience functions for creating synthetic coupled physiological signals
# useful for tutorials, vignettes, and quick testing.

#' Generate coupled sinusoidal signals
#'
#' Creates a pair of sinusoidal signals with a shared oscillatory component
#' at a specified coupling frequency and controllable noise level. Useful for
#' demonstrating and testing spectral coherence, phase synchrony, and other
#' cross-modal coupling measures.
#'
#' Both signals share the same sinusoidal component at \code{coupling_freq},
#' scaled by \code{coupling_strength}, with additive Gaussian noise scaled by
#' \code{noise}.
#'
#' @param sr1 Numeric sampling rate in Hz for the first signal (default 500).
#' @param sr2 Numeric sampling rate in Hz for the second signal (default 500).
#' @param coupling_freq Numeric frequency in Hz of the shared oscillatory
#'   component (default 20).
#' @param coupling_strength Numeric amplitude of the shared component
#'   (default 0.8).
#' @param noise Numeric standard deviation of the additive Gaussian noise
#'   (default 0.2).
#' @param duration Numeric duration in seconds (default 10).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{x}{Numeric vector -- first coupled signal.}
#'     \item{y}{Numeric vector -- second coupled signal.}
#'     \item{sr_x}{Numeric sampling rate of signal x.}
#'     \item{sr_y}{Numeric sampling rate of signal y.}
#'     \item{coupling_freq}{Numeric coupling frequency used.}
#'   }
#'
#' @references
#' Carter, G. C. (1987). Coherence and time delay estimation.
#' \emph{Proceedings of the IEEE}, 75(2), 236--255.
#'
#' @seealso [make_eeg_emg()], [make_directed_signals()],
#'   [coherence()], [couplingAnalysis()]
#' @export
#' @examples
#' signals <- make_coupled_signals(sr1 = 500, sr2 = 500,
#'                                  coupling_freq = 20,
#'                                  coupling_strength = 0.8,
#'                                  noise = 0.2, duration = 10)
#' result <- coherence(signals$x, signals$y, sr = signals$sr_x)
make_coupled_signals <- function(sr1 = 500, sr2 = 500,
                                  coupling_freq = 20,
                                  coupling_strength = 0.8,
                                  noise = 0.2, duration = 10) {

  stopifnot(is.numeric(sr1), sr1 > 0)
  stopifnot(is.numeric(sr2), sr2 > 0)
  stopifnot(is.numeric(coupling_freq), coupling_freq > 0)
  stopifnot(is.numeric(coupling_strength))
  stopifnot(is.numeric(noise), noise >= 0)
  stopifnot(is.numeric(duration), duration > 0)

  n1 <- as.integer(sr1 * duration)
  n2 <- as.integer(sr2 * duration)
  t1 <- seq(0, duration, length.out = n1)
  t2 <- seq(0, duration, length.out = n2)

  driver1 <- sin(2 * pi * coupling_freq * t1)
  driver2 <- sin(2 * pi * coupling_freq * t2)

  x <- coupling_strength * driver1 + noise * stats::rnorm(n1)
  y <- coupling_strength * driver2 + noise * stats::rnorm(n2)

  list(
    x = x,
    y = y,
    sr_x = sr1,
    sr_y = sr2,
    coupling_freq = coupling_freq
  )
}


#' Generate simulated EEG-EMG PhysioExperiment pair
#'
#' Creates a pair of \code{PhysioExperiment} objects simulating simultaneous
#' EEG and EMG recordings with corticomuscular coherence (CMC) at a specified
#' frequency. The first channel of each modality contains a shared oscillatory
#' component; remaining channels contain independent noise.
#'
#' This function is useful for demonstrating and testing multi-channel and
#' multi-modal coupling analyses such as coherence matrices and
#' \code{MultiPhysioExperiment} workflows.
#'
#' @param n_sec Numeric recording duration in seconds (default 10).
#' @param n_eeg_ch Integer number of EEG channels (default 3).
#' @param n_emg_ch Integer number of EMG channels (default 2).
#' @param sr_eeg Numeric EEG sampling rate in Hz (default 500).
#' @param sr_emg Numeric EMG sampling rate in Hz (default 1000).
#' @param cmc_freq Numeric corticomuscular coherence frequency in Hz
#'   (default 20).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{eeg}{A \code{PhysioExperiment} object with EEG data
#'       (\code{n_eeg_ch} channels at \code{sr_eeg} Hz).}
#'     \item{emg}{A \code{PhysioExperiment} object with EMG data
#'       (\code{n_emg_ch} channels at \code{sr_emg} Hz).}
#'   }
#'
#' @references
#' Halliday, D. M., Rosenberg, J. R., Amjad, A. M., Breeze, P.,
#' Conway, B. A., & Farmer, S. F. (1995). A framework for the analysis of
#' mixed time series/point process data -- theory and application to the
#' study of physiological tremor, single motor unit discharges and
#' electromyograms. \emph{Progress in Biophysics and Molecular Biology},
#' 64(2--3), 237--278.
#'
#' @seealso [make_coupled_signals()], [make_directed_signals()],
#'   [MultiPhysioExperiment()], [coherenceMatrix()]
#' @export
#' @examples
#' data <- make_eeg_emg(n_sec = 5, n_eeg_ch = 3, n_emg_ch = 2)
#' mpe <- MultiPhysioExperiment(EEG = data$eeg, EMG = data$emg)
#' modalities(mpe)
make_eeg_emg <- function(n_sec = 10, n_eeg_ch = 3, n_emg_ch = 2,
                          sr_eeg = 500, sr_emg = 1000, cmc_freq = 20) {

  stopifnot(is.numeric(n_sec), n_sec > 0)
  stopifnot(is.numeric(n_eeg_ch), n_eeg_ch >= 1)
  stopifnot(is.numeric(n_emg_ch), n_emg_ch >= 1)
  stopifnot(is.numeric(sr_eeg), sr_eeg > 0)
  stopifnot(is.numeric(sr_emg), sr_emg > 0)
  stopifnot(is.numeric(cmc_freq), cmc_freq > 0)

  n_eeg_ch <- as.integer(n_eeg_ch)
  n_emg_ch <- as.integer(n_emg_ch)

  n_eeg <- as.integer(sr_eeg * n_sec)
  n_emg <- as.integer(sr_emg * n_sec)
  t_eeg <- seq(0, n_sec, length.out = n_eeg)
  t_emg <- seq(0, n_sec, length.out = n_emg)

  driver_eeg <- sin(2 * pi * cmc_freq * t_eeg)
  driver_emg <- sin(2 * pi * cmc_freq * t_emg)

  eeg_data <- matrix(stats::rnorm(n_eeg * n_eeg_ch),
                     nrow = n_eeg, ncol = n_eeg_ch)
  eeg_data[, 1] <- eeg_data[, 1] + 0.7 * driver_eeg

  emg_data <- matrix(stats::rnorm(n_emg * n_emg_ch) * 0.5,
                     nrow = n_emg, ncol = n_emg_ch)
  emg_data[, 1] <- emg_data[, 1] + 0.7 * driver_emg

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


#' Generate directed coupling signals
#'
#' Creates a pair of signals where \code{x} drives \code{y} with a specified
#' lag and coupling strength. Signal \code{x} is white noise, and \code{y}
#' is a mixture of a lagged copy of \code{x} and independent noise. This is
#' useful for testing directed coupling measures such as Granger causality.
#'
#' The generating model is:
#' \deqn{y(t) = \text{coupling} \cdot x(t - \text{lag\_samples}) + (1 - \text{coupling}) \cdot \epsilon(t)}
#' where \eqn{\epsilon(t) \sim N(0, 1)}.
#'
#' @param n Integer number of samples (default 5000).
#' @param sr Numeric sampling rate in Hz (default 500).
#' @param lag_samples Integer number of samples by which \code{x} leads
#'   \code{y} (default 10).
#' @param coupling Numeric coupling strength in \[0, 1\] (default 0.7).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{x}{Numeric vector -- driving signal (white noise).}
#'     \item{y}{Numeric vector -- driven signal (lagged mixture).}
#'     \item{sr}{Numeric sampling rate.}
#'     \item{lag_samples}{Integer lag used.}
#'   }
#'
#' @references
#' Granger, C. W. J. (1969). Investigating causal relations by econometric
#' models and cross-spectral methods. \emph{Econometrica}, 37(3), 424--438.
#'
#' @seealso [grangerCausality()], [make_coupled_signals()],
#'   [crossCorrelation()], [couplingAnalysis()]
#' @export
#' @examples
#' signals <- make_directed_signals(n = 5000, sr = 500,
#'                                   lag_samples = 10, coupling = 0.7)
#' result <- grangerCausality(signals$x, signals$y, sr = signals$sr,
#'                            order = 15)
#' result$gc_xy   # should be positive (x drives y)
#' result$net_gc  # should be positive
make_directed_signals <- function(n = 5000, sr = 500,
                                   lag_samples = 10, coupling = 0.7) {

  stopifnot(is.numeric(n), n > 0)
  stopifnot(is.numeric(sr), sr > 0)
  stopifnot(is.numeric(lag_samples), lag_samples >= 1)
  stopifnot(is.numeric(coupling), coupling >= 0, coupling <= 1)

  n <- as.integer(n)
  lag_samples <- as.integer(lag_samples)

  x <- stats::rnorm(n)
  y <- numeric(n)
  for (i in (lag_samples + 1L):n) {
    y[i] <- coupling * x[i - lag_samples] + (1 - coupling) * stats::rnorm(1)
  }

  list(x = x, y = y, sr = sr, lag_samples = lag_samples)
}
