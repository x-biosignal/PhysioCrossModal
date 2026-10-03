# PhysioCrossModal 0.8.4

## Documentation

* `?PhysioCrossModal` now answers: a package help page gives one paragraph on what the
  package is for, the main entry points grouped by task, and where to go next.
* Runnable `@examples` added or corrected across 1 help pages. Each runs
  offline in seconds, writes nothing outside `tempdir()`, and is executed by
  `R CMD check`; anything needing a device, a download or an optional backend is
  fenced with the reason stated.
* The README's quick start runs as written: it attaches the package, builds its
  own inputs, and uses only hard dependencies.

# PhysioCrossModal 0.8.3

- `.reject_gap_in_window()` no longer carries its own gap detection. It calls
  `PhysioExperiment::detectGaps()`, the same function `streamCoverage()` judges
  trial completeness with. Sharing the 1.5x factor was not enough to make the two
  agree: this guard detected on the whole series and then asked which gaps met the
  window, while the foundation selected the window's samples first and
  differenced those, which loses a gap straddling an edge. For a stream with a
  0.2 s hole at 10 Hz, a window of [0.25, 0.45] was refused here and reported as
  fully covered there. One implementation, one answer. Requires
  PhysioExperiment >= 2.1.3.

# PhysioCrossModal 0.8.2

* **The pinned Tensorpac fixture now ships with the tests.** `inst/validation` is
  stripped from the published package by design, so 0.8.1 shipped this test
  without the data it reads and `R CMD check` failed on **seven** builds -- Linux
  release and devel, Windows release, devel and oldrel, and macOS release and
  oldrel. The source package and the Wasm build succeeded, which is why 0.8.1
  reached the distribution index with `check: ERROR` recorded against its
  binaries.

  The fixture moved to `tests/testthat/fixtures/` (21 kB: input values, expected
  values, checksums, provenance). It contains no Tensorpac code and runs no
  Python -- the test compares this package's computation against pinned reference
  numbers -- so agreement with the external reference is now verified in the
  published package rather than only in the development tree. The generator that
  produced the numbers stays unpublished in `inst/validation`, and `SHA256SUMS` is
  scoped to the files that ship, with the generator's identity still pinned by
  `Generator-SHA256` in the manifest.

# PhysioCrossModal 0.8.1

## Real sample times are consulted, and the cache is keyed on them

* **An equal start offset is no longer treated as an equal time
  correspondence.** `.extract_signal_pair()` only reached for real sample times
  when the offsets differed, so two modalities that began together but diverged
  afterwards -- a gap, a drifting clock -- were still paired element for element.
  With a one-second hole in one stream, t = 3.57 s was compared against
  t = 2.57 s. The historical path is now taken only when both streams are on a
  reconstructed regular grid from a common start, which is the case whose
  existing numbers must not move.

* **Interpolating across a gap is refused rather than smoothed over.** A step
  between measured sample times larger than 1.5x the nominal period is treated
  as a gap; a gap inside the requested window is an error naming the gap and how
  many samples it spans. 1.5 is an adopted criterion that separates acquisition
  jitter from a dropped sample -- it is not a proof of absence under every
  condition, and the boundaries are documented at `.reject_gap_in_window()`.
  A gap outside the window is no obstacle.

* **The coupling cache is keyed on the measured times.** The fingerprint covered
  signal values, rate, `t0` and offset but not the times that decide which
  moments are compared, so changing only the times returned a stale result --
  differing by as much as 0.91 in coherence. Changing them now misses the cache
  and recomputes; identical input still hits it.

* `coherence()` reports the alignment its estimate rests on under `$input`, and
  the fast-path flag is named `regular_common_start` rather than the earlier
  `identical_grids`, which overclaimed: at 100 Hz beside 200 Hz the time vectors
  are not identical.

# PhysioCrossModal 0.8.0

## Breaking: the container is the ecosystem's, not this package's

* This package no longer defines its own `MultiPhysioExperiment`. The one
  multimodal container lives in `PhysioExperiment` and is re-exported here, so
  the constructor and every accessor (`experiments()`, `modalities()`,
  `samplingRates()`, `nModalities()`, `alignment()`, `alignment<-`) keep their
  names and return types. Two S4 classes used to share this name, and whichever
  package loaded second silently redefined the other.

* **Coupling between modalities that did not start together is now computed on
  their overlap.** `.extract_signal_pair()` matched sampling rates but then
  compared the two arrays head to head, ignoring the start offsets the container
  records, so modalities recorded a second apart were estimated over different
  time spans. Results on simultaneously started data are unchanged; a test
  asserts this against the raw assay columns. Modalities that do not overlap at
  all are now an error naming both spans.

* **The coupling cache is keyed on its inputs.** `couplingAnalysisCached()` keyed
  only on the method, modality names, channel indices and extra arguments, so
  editing the signals or shifting the clock and re-running returned the previous
  answer. The key now covers the signal values, assay, rate, clock, channels,
  method, parameters and implementation version. The cache is session-scoped
  rather than stored in the container, so it is no longer carried inside saved
  objects; recompute after loading.

* An empty container is now allowed, `x[c(a, a), ]` selects a one-sample window
  rather than erroring, and error messages speak of *streams*.

## Unreleased notes (pre-0.7.1)

- New `alignObservationBlocks()` validates explicit composite observation keys, aligns feature blocks, and reports row mappings and any explicitly requested intersection exclusions before positional fusion.

# PhysioCrossModal 0.7.1

- `mutualInformation()` — mutual information (Shannon 1948) between two signals via the histogram
  estimator (equal-width bins, in nats): the symmetric, undirected, nonlinear measure of statistical
  dependence, the information-theoretic complement of `transferEntropy()` (directed). Reproduces
  `sklearn.metrics.mutual_info_score` bit-for-bit on the same binning. The histogram estimator is
  positively biased, so compare against a shuffled control.

# PhysioCrossModal 0.6.4

- `leidaTransitions()` — LEiDA state-transition metrics: the Markov transition
  matrix between connectivity states, switching rate, fractional occupancy, mean
  dwell time, and occupancy-weighted transition entropy.

# PhysioCrossModal 0.6.3

- `crossWaveletTransform()` — standalone cross-wavelet transform (XWT) with
  time-frequency common power, phase lead/lag, cone of influence, and AR(1)
  red-noise significance by Monte-Carlo surrogates (Torrence & Compo 1998;
  Grinsted et al. 2004).

# PhysioCrossModal 0.6.2

- New frontier connectivity measures (`R/coupling-frontier.R`):
  - `phaseSlopeIndex()` — directed, volume-conduction-robust flow direction
    (Nolte et al. 2008); positive means the first signal leads.
  - `orthogonalizedAEC()` — leakage-corrected amplitude-envelope correlation
    (Hipp et al. 2012; Brookes et al. 2012).
  - `ciPLV()` — corrected imaginary phase-locking value (Bruna et al. 2018).
  - `pairwisePhaseConsistency()` — bias-free phase synchrony (Vinck et al. 2010).
  - `leidaStates()` — LEiDA dynamic functional-connectivity states
    (Cabral et al. 2017).

# PhysioCrossModal 0.6.1

- Fixed `elasticAlign()` / `srvfMean()` under `fdasrvf` >= 2.x: the aligned
  function is now read from `fit$f2tilde` (renamed from `f2n` upstream), with a
  fallback to `f2n` for older `fdasrvf`. Previously the aligned curve came back
  empty on newer `fdasrvf`, breaking the elastic mean's length invariant (this
  surfaced only when `fdasrvf` was installed, e.g. on r-universe binary builds).
- Replaced an unexercised `elasticAlign()` test assertion (`amplitude_distance <
  1e-3` on the `fdasrvf` path, which was skipped whenever `fdasrvf` was absent)
  with realistic checks: the aligned curve tracks the reference and the
  amplitude distance is reduced by more than an order of magnitude.

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
