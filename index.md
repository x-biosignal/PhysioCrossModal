# PhysioCrossModal ![PhysioCrossModal logo](reference/figures/logo.png)

**Cross-Modal Coupling Analysis for PhysioExperiment Objects**

PhysioCrossModal provides cross-modal coupling, connectivity, and
synchrony analysis between physiological signals of different modalities
(EEG, EMG, ECG, EDA, MoCap, fNIRS, etc.). It introduces the
`MultiPhysioExperiment` container class for holding multiple
`PhysioExperiment` objects at different sampling rates with temporal
alignment. Its exported analysis surface supports spectral coherence,
phase synchrony (PLV, PLI, wPLI), phase-amplitude coupling, directed
coupling (Granger causality), time-domain coupling (cross-correlation),
surrogate-based statistical testing, and publication-ready
visualization.

## Installation

You can install PhysioCrossModal from
[r-universe](https://x-biosignal.r-universe.dev):

``` r

install.packages("PhysioCrossModal",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

Or install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioCrossModal")
```

## Quick Start

``` r

library(PhysioCrossModal)

# Create simulated EEG + EMG data with known corticomuscular coupling
sim <- make_eeg_emg(
  n_time = 5000,
  eeg_rate = 500,
  emg_rate = 1000,
  coupling_freq = 20   # beta-band coupling
)

# Build a MultiPhysioExperiment from the two modalities
mpe <- MultiPhysioExperiment(
  experiments = list(EEG = sim$eeg, EMG = sim$emg)
)

# Inspect the container
mpe
modalities(mpe)      # c("EEG", "EMG")
nModalities(mpe)     # 2
samplingRates(mpe)   # named vector of sampling rates

# Align signals to a common time base
mpe_aligned <- alignSignals(mpe)

# Compute corticomuscular coherence
coh <- coherence(
  mpe_aligned,
  x_modality = "EEG",
  y_modality = "EMG",
  method = "multitaper"
)

# Visualize the coherence spectrum
plotCoherenceSpectrum(coh, freq_range = c(1, 50))

# Test statistical significance with phase-shuffled surrogates
sig <- surrogateTest(coh, n_surrogates = 200)
sig$p_value
```

## Features

### MultiPhysioExperiment Class

The `MultiPhysioExperiment` S4 class is a container for multiple
`PhysioExperiment` objects recorded at different sampling rates. It
maintains temporal alignment metadata so that cross-modal analyses
operate on correctly synchronized signals.

- **Constructor:** `MultiPhysioExperiment(experiments, alignment)`
- **Accessors:**
  [`experiments()`](https://x-biosignal.github.io/PhysioCrossModal/reference/experiments.md),
  [`alignment()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignment.md),
  [`modalities()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modalities.md),
  [`nModalities()`](https://x-biosignal.github.io/PhysioCrossModal/reference/nModalities.md),
  [`samplingRates()`](https://x-biosignal.github.io/PhysioCrossModal/reference/samplingRates.md)
- **Manipulation:**
  [`mergePhysio()`](https://x-biosignal.github.io/PhysioCrossModal/reference/mergePhysio.md)
  – combine two MultiPhysioExperiment objects
- **Subsetting:** `[` and `[[` methods for extracting modalities by name
  or index

### Signal Alignment

Temporal synchronization across modalities with different sampling
rates:

- [`alignSignals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignSignals.md)
  – align all modalities to a common time base using interpolation
- [`alignToRate()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignToRate.md)
  – resample all modalities to a specified target sampling rate

### Spectral Coupling

Frequency-domain coupling measures for identifying shared oscillatory
activity:

- [`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md)
  – magnitude-squared coherence between signal pairs
- [`coherenceMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherenceMatrix.md)
  – pairwise coherence across all channel combinations
- [`crossSpectrum()`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossSpectrum.md)
  – complex-valued cross-spectral density
- [`multitaperCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/multitaperCoherence.md)
  – coherence estimated via multitaper spectral analysis (DPSS tapers)
- [`waveletCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletCoherence.md)
  – time-frequency coherence using continuous wavelet transform

### Phase Synchrony

Phase-based coupling measures robust to amplitude variations:

- [`phaseLockingValue()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLockingValue.md)
  – PLV for quantifying inter-trial or inter-signal phase consistency
- [`phaseLagIndex()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLagIndex.md)
  – PLI, robust to volume conduction artifacts
- [`weightedPLI()`](https://x-biosignal.github.io/PhysioCrossModal/reference/weightedPLI.md)
  – wPLI, weighted variant with improved sensitivity
- [`waveletPLV()`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletPLV.md)
  – time-resolved PLV using wavelet-extracted instantaneous phase

### Phase-Amplitude Coupling

- [`phaseAmplitudeCoupling()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseAmplitudeCoupling.md)
  – one declared phase/amplitude band using Tort, Canolty, Ozkurt, or
  PLV statistics
- [`comodulogram()`](https://x-biosignal.github.io/PhysioCrossModal/reference/comodulogram.md)
  – observed PAC across a phase-frequency by amplitude-frequency grid
- [`modulationIndex()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modulationIndex.md)
  – strict Tort/Canolty/Ozkurt comodulogram with one surrogate signal
  per complete grid, conservative cell-wise p-values, and grid-wide
  multiplicity adjustment

PAC is an association, not evidence of a biological mechanism. Waveform
shape, filter leakage and transients, non-stationarity, and common
inputs can all produce high or significant values. Inspect the realized
pass bands and the returned diagnostic strings before interpreting a
significance mask.

### Directed Coupling

Methods for inferring directional (causal) interactions between signals:

- [`grangerCausality()`](https://x-biosignal.github.io/PhysioCrossModal/reference/grangerCausality.md)
  – parametric time/frequency-domain Granger causality plus optional
  nonparametric Welch/DPSS estimation with Wilson spectral factorization
  and frequency-resolved Geweke measures

Granger causality is a stationary linear predictive-direction measure.
It does not establish physical or interventional causality and does not
remove confounding, common-driver, volume-conduction, or
instantaneous-mixing bias.

### Time-Domain Coupling

Coupling measures operating directly on the time series:

- [`crossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossCorrelation.md)
  – normalized cross-correlation with optimal lag estimation
- [`slidingCrossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/slidingCrossCorrelation.md)
  – time-resolved cross-correlation in sliding windows

### Unified Coupling API

High-level wrappers that dispatch to specific coupling methods:

- [`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)
  – compute any coupling measure via a single interface with `method`
  argument
- [`couplingMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingMatrix.md)
  – pairwise coupling matrix across all channel pairs for any supported
  method

### Statistical Inference

Non-parametric significance testing and generalization assessment:

- [`surrogateTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateTest.md)
  – phase-shuffled surrogate testing for a single coupling value
- [`surrogateMatrixTest()`](https://x-biosignal.github.io/PhysioCrossModal/reference/surrogateMatrixTest.md)
  – surrogate testing for full coupling matrices with multiple
  comparison correction
- [`modulationIndex()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modulationIndex.md)
  – full-grid PAC surrogates with adjusted cell-wise significance
- [`bootstrapCI()`](https://x-biosignal.github.io/PhysioCrossModal/reference/bootstrapCI.md)
  – bootstrap confidence intervals for coupling estimates
- [`lodoGeneralization()`](https://x-biosignal.github.io/PhysioCrossModal/reference/lodoGeneralization.md)
  – leave-one-dataset-out cross-validation for assessing
  generalizability of coupling findings across datasets or subjects

### Visualization

Publication-ready plotting functions for coupling results:

- [`plotCoherenceSpectrum()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotCoherenceSpectrum.md)
  – coherence as a function of frequency with significance threshold
- [`plotCouplingMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotCouplingMatrix.md)
  – heatmap of pairwise coupling values across channels
- [`plotCouplingTimecourse()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotCouplingTimecourse.md)
  – time-resolved coupling (e.g., sliding cross-correlation or wavelet
  PLV)
- [`plotWaveletCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotWaveletCoherence.md)
  – time-frequency coherence scalogram
- [`plotComodulogram()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotComodulogram.md)
  – PAC heatmap with optional significance attenuation

### Simulated Data

Helper functions for generating test data with known coupling
properties:

- [`make_coupled_signals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_coupled_signals.md)
  – two signals with controlled spectral coupling
- [`make_directed_signals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_directed_signals.md)
  – signals with known causal (directed) relationship
- [`make_eeg_emg()`](https://x-biosignal.github.io/PhysioCrossModal/reference/make_eeg_emg.md)
  – realistic EEG and EMG pair with corticomuscular coherence

## Use Cases

| Application | Methods |
|----|----|
| Corticomuscular coherence (EEG-EMG) | [`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md), [`multitaperCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/multitaperCoherence.md) |
| EEG inter-regional connectivity | [`phaseLockingValue()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLockingValue.md), [`phaseLagIndex()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseLagIndex.md), [`weightedPLI()`](https://x-biosignal.github.io/PhysioCrossModal/reference/weightedPLI.md) |
| Brain-heart coupling (EEG-ECG) | [`crossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossCorrelation.md), [`waveletCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletCoherence.md) |
| Neurovascular coupling (EEG-fNIRS) | [`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md), [`grangerCausality()`](https://x-biosignal.github.io/PhysioCrossModal/reference/grangerCausality.md) |
| EDA-EMG sympathetic co-activation | [`slidingCrossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/slidingCrossCorrelation.md), [`couplingMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingMatrix.md) |
| Movement-brain synchrony (MoCap-EEG) | [`waveletPLV()`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletPLV.md), [`waveletCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/waveletCoherence.md) |

## Dependencies

- **R** (\>= 4.2)
- **PhysioCore** – core data structures and accessors
- **methods**, **stats** – base R infrastructure
- **SummarizedExperiment**, **S4Vectors** – Bioconductor infrastructure

Optional (in Suggests):

- **signal** – additional DSP functions
- **PhysioPreprocess** – signal preprocessing
- **ggplot2** – visualization
- **knitr**, **rmarkdown** – vignettes

## PhysioExperiment Ecosystem

PhysioCrossModal is part of the PhysioExperiment ecosystem, a suite of R
packages for multi-modal physiological signal analysis:

| Package | Description |
|----|----|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| [PhysioPreprocess](https://github.com/x-biosignal/PhysioPreprocess) | Preprocessing (filters, ICA, resampling) |
| [PhysioAnalysis](https://github.com/x-biosignal/PhysioAnalysis) | Analysis and visualization |
| **PhysioCrossModal** | Cross-modal coupling and connectivity |
| [PhysioMoCap](https://github.com/x-biosignal/PhysioMoCap) | Motion capture and biomechanics |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to
browse all available packages.

## License

MIT License. See
[LICENSE](https://x-biosignal.github.io/PhysioCrossModal/LICENSE) for
details.

## Author

Yusuke Matsui

## Governance & support

Part of the [Physio ecosystem](https://x-biosignal.r-universe.dev).
Community and policy documents live in the umbrella repository:

- [Code of
  Conduct](https://github.com/x-biosignal/PhysioExperiment/blob/main/CODE_OF_CONDUCT.md)
- [Contributing](https://github.com/x-biosignal/PhysioExperiment/blob/main/CONTRIBUTING.md)
- [Governance](https://github.com/x-biosignal/PhysioExperiment/blob/main/GOVERNANCE.md)
- [Support](https://github.com/x-biosignal/PhysioExperiment/blob/main/SUPPORT.md)
- [Security
  policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/SECURITY.md)
- [Deprecation & lifecycle
  policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/DEPRECATION.md)
