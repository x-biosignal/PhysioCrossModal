# Changelog

## PhysioCrossModal 0.6.1

- Fixed
  [`elasticAlign()`](https://x-biosignal.github.io/PhysioCrossModal/reference/elasticAlign.md)
  /
  [`srvfMean()`](https://x-biosignal.github.io/PhysioCrossModal/reference/srvfMean.md)
  under `fdasrvf` \>= 2.x: the aligned function is now read from
  `fit$f2tilde` (renamed from `f2n` upstream), with a fallback to `f2n`
  for older `fdasrvf`. Previously the aligned curve came back empty on
  newer `fdasrvf`, breaking the elastic mean’s length invariant (this
  surfaced only when `fdasrvf` was installed, e.g. on r-universe binary
  builds).
- Replaced an unexercised
  [`elasticAlign()`](https://x-biosignal.github.io/PhysioCrossModal/reference/elasticAlign.md)
  test assertion (`amplitude_distance < 1e-3` on the `fdasrvf` path,
  which was skipped whenever `fdasrvf` was absent) with realistic
  checks: the aligned curve tracks the reference and the amplitude
  distance is reduced by more than an order of magnitude.

## PhysioCrossModal 0.6.0

- Added
  [`modulationIndex()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modulationIndex.md)
  as the inference-enabled compatibility API for Tort, Canolty, and
  Ozkurt phase-amplitude comodulograms. It validates exact pass bands,
  reuses filtered phase/amplitude components, and returns cell-wise
  conservative surrogate p-values, grid-wide multiplicity adjustment,
  thresholds, and an auditable significance mask.
- Surrogate replicates are generated once for the complete frequency
  grid before deterministic sequential or parallel evaluation. Seeded
  results are core-independent and the caller’s random-number state is
  preserved.
- [`plotComodulogram()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotComodulogram.md)
  now accepts validated logical or p-value masks and attenuates
  non-significant cells without changing observed PAC values.
- PAC documentation now distinguishes association from mechanism and
  calls out waveform, filter, edge, non-stationarity, and common-input
  artifacts.

## PhysioCrossModal 0.5.0

- Added nonparametric bivariate Granger causality using Welch or DPSS
  multitaper cross-spectral estimation, Wilson spectral factorization,
  and Geweke directional spectra. Results include convergence,
  reconstruction, regularization, and estimator diagnostics.
- Added the named-only `granger_method` selector to
  [`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)
  so the unified dispatcher’s coupling-family argument no longer
  collides with the Granger estimator selector.
- The existing parametric Granger implementation remains the default and
  retains its numerical behavior.

## PhysioCrossModal 0.4.0

Initial release as a standalone package in the x-biosignal /
PhysioExperiment ecosystem. PhysioCrossModal provides cross-modal
coupling, connectivity, and synchrony analysis between physiological
signals of different modalities (EEG, EMG, ECG, EDA, MoCap, fNIRS,
etc.).

### New Features

- Added the `MultiPhysioExperiment` S4 container for holding several
  `PhysioExperiment` objects recorded simultaneously at potentially
  different sampling rates, with temporal alignment metadata and a
  coupling-result cache:
  - Accessors
    [`experiments()`](https://x-biosignal.github.io/PhysioCrossModal/reference/experiments.md),
    [`modalities()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modalities.md),
    [`nModalities()`](https://x-biosignal.github.io/PhysioCrossModal/reference/nModalities.md),
    [`samplingRates()`](https://x-biosignal.github.io/PhysioCrossModal/reference/samplingRates.md),
    [`alignment()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignment.md),
    and
    [`couplingResults()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingResults.md)
    (with replacement forms) plus `[[`, `length`, `names`, and a `show`
    method.
- Added signal alignment and merging utilities:
  - [`alignToRate()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignToRate.md)
    resamples a `PhysioExperiment` to a target rate via linear, spline,
    or FFT interpolation, with anti-alias lowpass filtering before
    downsampling.
  - [`alignSignals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignSignals.md)
    brings several modalities onto a common rate (`lowest_rate`,
    `common_rate`, or explicit `resample`) and wraps them in a
    `MultiPhysioExperiment`.
  - [`mergePhysio()`](https://x-biosignal.github.io/PhysioCrossModal/reference/mergePhysio.md)
    concatenates channels of two rate-matched objects with label
    prefixing.
- Added spectral coupling:
  [`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md)
  and
  [`multitaperCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/multitaperCoherence.md)
  (magnitude-squared coherence via Welch and Thomson multitaper
  estimation) and
  [`crossSpectrum()`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossSpectrum.md)
  (cross-spectral density).
- Added phase-synchrony measures over a chosen frequency band:
  [`phaseLockingValue()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLockingValue.md)
  (PLV),
  [`phaseLagIndex()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLagIndex.md)
  (PLI), and
  [`weightedPLI()`](https://x-biosignal.github.io/PhysioCrossModal/reference/weightedPLI.md)
  (wPLI, with optional debiasing).
- Added directed coupling via
  [`grangerCausality()`](https://x-biosignal.github.io/PhysioCrossModal/reference/grangerCausality.md),
  supporting both time-domain (parametric VAR) and Geweke spectral
  Granger causality.
- Added time-domain coupling:
  [`crossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossCorrelation.md)
  with lag estimation and
  [`slidingCrossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/slidingCrossCorrelation.md)
  for tracking time-varying coupling.
- Added time-frequency coupling with complex Morlet wavelets:
  [`waveletCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletCoherence.md)
  and
  [`waveletPLV()`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletPLV.md),
  including a cone-of-influence estimate.
- Added a unified dispatcher
  [`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)
  that routes numeric vectors, `PhysioExperiment` pairs, or
  `MultiPhysioExperiment` modalities to any coupling method, plus
  [`couplingAnalysisCached()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysisCached.md)
  for `digest`-keyed memoisation into the object’s `couplingResults`
  cache.
- Added multi-channel coupling matrices
  [`coherenceMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherenceMatrix.md)
  and
  [`couplingMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingMatrix.md)
  that compute a statistic for every channel pair.

### Statistical Significance

- Added surrogate-based significance testing with
  [`surrogateTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateTest.md),
  using phase-randomisation or time-shift surrogates and the
  conservative Phipson & Smyth p-value correction, with optional
  multi-core execution.
- Added
  [`surrogateMatrixTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateMatrixTest.md)
  for pairwise significance across a coupling matrix with FDR
  correction, and
  [`bootstrapCI()`](https://x-biosignal.github.io/PhysioCrossModal/reference/bootstrapCI.md)
  for moving-block bootstrap confidence intervals of coupling
  statistics.
- Added
  [`lodoGeneralization()`](https://x-biosignal.github.io/PhysioCrossModal/reference/lodoGeneralization.md),
  a leave-one-site-out benchmark for evaluating model transportability
  across sites for Gaussian and binomial outcomes.

### Visualization

- Added `ggplot2`-based plots:
  [`plotCouplingMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotCouplingMatrix.md)
  (heatmap),
  [`plotCoherenceSpectrum()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotCoherenceSpectrum.md),
  [`plotCouplingTimecourse()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotCouplingTimecourse.md)
  (sliding-window), and
  [`plotWaveletCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotWaveletCoherence.md)
  (time-frequency map).

### Utilities

- Added simulated-data generators for testing and demonstration:
  [`make_coupled_signals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_coupled_signals.md),
  [`make_eeg_emg()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_eeg_emg.md)
  (corticomuscular coherence), and
  [`make_directed_signals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_directed_signals.md).
