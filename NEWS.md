# PhysioCrossModal 0.6.0

- Added `modulationIndex()` as the inference-enabled compatibility API for
  Tort, Canolty, and Ozkurt phase-amplitude comodulograms. It validates exact
  pass bands, reuses filtered phase/amplitude components, and returns
  cell-wise conservative surrogate p-values, grid-wide multiplicity
  adjustment, thresholds, and an auditable significance mask.
- Surrogate replicates are generated once for the complete frequency grid
  before deterministic sequential or parallel evaluation. Seeded results are
  core-independent and the caller's random-number state is preserved.
- `plotComodulogram()` now accepts validated logical or p-value masks and
  attenuates non-significant cells without changing observed PAC values.
- PAC documentation now distinguishes association from mechanism and calls out
  waveform, filter, edge, non-stationarity, and common-input artifacts.

# PhysioCrossModal 0.5.0

- Added nonparametric bivariate Granger causality using Welch or DPSS
  multitaper cross-spectral estimation, Wilson spectral factorization, and
  Geweke directional spectra. Results include convergence, reconstruction,
  regularization, and estimator diagnostics.
- Added the named-only `granger_method` selector to `couplingAnalysis()` so
  the unified dispatcher's coupling-family argument no longer collides with
  the Granger estimator selector.
- The existing parametric Granger implementation remains the default and
  retains its numerical behavior.

# PhysioCrossModal 0.4.0

Initial release as a standalone package in the x-biosignal / PhysioExperiment
ecosystem. PhysioCrossModal provides cross-modal coupling, connectivity, and
synchrony analysis between physiological signals of different modalities
(EEG, EMG, ECG, EDA, MoCap, fNIRS, etc.).

## New Features

- Added the `MultiPhysioExperiment` S4 container for holding several
  `PhysioExperiment` objects recorded simultaneously at potentially different
  sampling rates, with temporal alignment metadata and a coupling-result cache:
  - Accessors `experiments()`, `modalities()`, `nModalities()`,
    `samplingRates()`, `alignment()`, and `couplingResults()` (with
    replacement forms) plus `[[`, `length`, `names`, and a `show` method.
- Added signal alignment and merging utilities:
  - `alignToRate()` resamples a `PhysioExperiment` to a target rate via
    linear, spline, or FFT interpolation, with anti-alias lowpass filtering
    before downsampling.
  - `alignSignals()` brings several modalities onto a common rate
    (`lowest_rate`, `common_rate`, or explicit `resample`) and wraps them in a
    `MultiPhysioExperiment`.
  - `mergePhysio()` concatenates channels of two rate-matched objects with
    label prefixing.
- Added spectral coupling: `coherence()` and `multitaperCoherence()`
  (magnitude-squared coherence via Welch and Thomson multitaper estimation)
  and `crossSpectrum()` (cross-spectral density).
- Added phase-synchrony measures over a chosen frequency band:
  `phaseLockingValue()` (PLV), `phaseLagIndex()` (PLI), and `weightedPLI()`
  (wPLI, with optional debiasing).
- Added directed coupling via `grangerCausality()`, supporting both
  time-domain (parametric VAR) and Geweke spectral Granger causality.
- Added time-domain coupling: `crossCorrelation()` with lag estimation and
  `slidingCrossCorrelation()` for tracking time-varying coupling.
- Added time-frequency coupling with complex Morlet wavelets:
  `waveletCoherence()` and `waveletPLV()`, including a cone-of-influence
  estimate.
- Added a unified dispatcher `couplingAnalysis()` that routes numeric vectors,
  `PhysioExperiment` pairs, or `MultiPhysioExperiment` modalities to any
  coupling method, plus `couplingAnalysisCached()` for `digest`-keyed
  memoisation into the object's `couplingResults` cache.
- Added multi-channel coupling matrices `coherenceMatrix()` and
  `couplingMatrix()` that compute a statistic for every channel pair.

## Statistical Significance

- Added surrogate-based significance testing with `surrogateTest()`, using
  phase-randomisation or time-shift surrogates and the conservative
  Phipson & Smyth p-value correction, with optional multi-core execution.
- Added `surrogateMatrixTest()` for pairwise significance across a coupling
  matrix with FDR correction, and `bootstrapCI()` for moving-block bootstrap
  confidence intervals of coupling statistics.
- Added `lodoGeneralization()`, a leave-one-site-out benchmark for evaluating
  model transportability across sites for Gaussian and binomial outcomes.

## Visualization

- Added `ggplot2`-based plots: `plotCouplingMatrix()` (heatmap),
  `plotCoherenceSpectrum()`, `plotCouplingTimecourse()` (sliding-window), and
  `plotWaveletCoherence()` (time-frequency map).

## Utilities

- Added simulated-data generators for testing and demonstration:
  `make_coupled_signals()`, `make_eeg_emg()` (corticomuscular coherence), and
  `make_directed_signals()`.
