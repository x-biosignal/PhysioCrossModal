# PhysioCrossModal <img src="man/figures/logo.png" align="right" height="139" alt="PhysioCrossModal logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/x-biosignal/PhysioCrossModal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/x-biosignal/PhysioCrossModal/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/PhysioCrossModal)](https://CRAN.R-project.org/package=PhysioCrossModal)
[![r-universe](https://x-biosignal.r-universe.dev/badges/PhysioCrossModal)](https://x-biosignal.r-universe.dev/PhysioCrossModal)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Cross-Modal Coupling Analysis for PhysioExperiment Objects**

PhysioCrossModal provides cross-modal coupling, connectivity, and synchrony analysis between physiological signals of different modalities (EEG, EMG, ECG, EDA, MoCap, fNIRS, etc.). It introduces the `MultiPhysioExperiment` container class for holding multiple `PhysioExperiment` objects at different sampling rates with temporal alignment. Its exported analysis surface supports spectral coherence, phase synchrony (PLV, PLI, wPLI), phase-amplitude coupling, directed coupling (Granger causality), time-domain coupling (cross-correlation), surrogate-based statistical testing, and publication-ready visualization.

## Installation

You can install PhysioCrossModal from [r-universe](https://x-biosignal.r-universe.dev):

```r
install.packages("PhysioCrossModal",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioCrossModal")
```

## Quick Start

```r
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

The `MultiPhysioExperiment` S4 class is a container for multiple `PhysioExperiment` objects recorded at different sampling rates. It maintains temporal alignment metadata so that cross-modal analyses operate on correctly synchronized signals.

- **Constructor:** `MultiPhysioExperiment(experiments, alignment)`
- **Accessors:** `experiments()`, `alignment()`, `modalities()`, `nModalities()`, `samplingRates()`
- **Manipulation:** `mergePhysio()` -- combine two MultiPhysioExperiment objects
- **Subsetting:** `[` and `[[` methods for extracting modalities by name or index

### Signal Alignment

Temporal synchronization across modalities with different sampling rates:

- `alignSignals()` -- align all modalities to a common time base using interpolation
- `alignToRate()` -- resample all modalities to a specified target sampling rate

### Spectral Coupling

Frequency-domain coupling measures for identifying shared oscillatory activity:

- `coherence()` -- magnitude-squared coherence between signal pairs
- `coherenceMatrix()` -- pairwise coherence across all channel combinations
- `crossSpectrum()` -- complex-valued cross-spectral density
- `multitaperCoherence()` -- coherence estimated via multitaper spectral analysis (DPSS tapers)
- `waveletCoherence()` -- time-frequency coherence using continuous wavelet transform

### Phase Synchrony

Phase-based coupling measures robust to amplitude variations:

- `phaseLockingValue()` -- PLV for quantifying inter-trial or inter-signal phase consistency
- `phaseLagIndex()` -- PLI, robust to volume conduction artifacts
- `weightedPLI()` -- wPLI, weighted variant with improved sensitivity
- `waveletPLV()` -- time-resolved PLV using wavelet-extracted instantaneous phase

### Phase-Amplitude Coupling

- `phaseAmplitudeCoupling()` -- one declared phase/amplitude band using Tort,
  Canolty, Ozkurt, or PLV statistics
- `comodulogram()` -- observed PAC across a phase-frequency by
  amplitude-frequency grid
- `modulationIndex()` -- strict Tort/Canolty/Ozkurt comodulogram with one
  surrogate signal per complete grid, conservative cell-wise p-values, and
  grid-wide multiplicity adjustment

PAC is an association, not evidence of a biological mechanism. Waveform shape,
filter leakage and transients, non-stationarity, and common inputs can all
produce high or significant values. Inspect the realized pass bands and the
returned diagnostic strings before interpreting a significance mask.

### Directed Coupling

Methods for inferring directional (causal) interactions between signals:

- `grangerCausality()` -- parametric time/frequency-domain Granger causality
  plus optional nonparametric Welch/DPSS estimation with Wilson spectral
  factorization and frequency-resolved Geweke measures

Granger causality is a stationary linear predictive-direction measure. It
does not establish physical or interventional causality and does not remove
confounding, common-driver, volume-conduction, or instantaneous-mixing bias.

### Time-Domain Coupling

Coupling measures operating directly on the time series:

- `crossCorrelation()` -- normalized cross-correlation with optimal lag estimation
- `slidingCrossCorrelation()` -- time-resolved cross-correlation in sliding windows

### Unified Coupling API

High-level wrappers that dispatch to specific coupling methods:

- `couplingAnalysis()` -- compute any coupling measure via a single interface with `method` argument
- `couplingMatrix()` -- pairwise coupling matrix across all channel pairs for any supported method

### Statistical Inference

Non-parametric significance testing and generalization assessment:

- `surrogateTest()` -- phase-shuffled surrogate testing for a single coupling value
- `surrogateMatrixTest()` -- surrogate testing for full coupling matrices with multiple comparison correction
- `modulationIndex()` -- full-grid PAC surrogates with adjusted cell-wise
  significance
- `bootstrapCI()` -- bootstrap confidence intervals for coupling estimates
- `lodoGeneralization()` -- leave-one-dataset-out cross-validation for assessing generalizability of coupling findings across datasets or subjects

### Visualization

Publication-ready plotting functions for coupling results:

- `plotCoherenceSpectrum()` -- coherence as a function of frequency with significance threshold
- `plotCouplingMatrix()` -- heatmap of pairwise coupling values across channels
- `plotCouplingTimecourse()` -- time-resolved coupling (e.g., sliding cross-correlation or wavelet PLV)
- `plotWaveletCoherence()` -- time-frequency coherence scalogram
- `plotComodulogram()` -- PAC heatmap with optional significance attenuation

### Simulated Data

Helper functions for generating test data with known coupling properties:

- `make_coupled_signals()` -- two signals with controlled spectral coupling
- `make_directed_signals()` -- signals with known causal (directed) relationship
- `make_eeg_emg()` -- realistic EEG and EMG pair with corticomuscular coherence

## Use Cases

| Application | Methods |
|-------------|---------|
| Corticomuscular coherence (EEG-EMG) | `coherence()`, `multitaperCoherence()` |
| EEG inter-regional connectivity | `phaseLockingValue()`, `phaseLagIndex()`, `weightedPLI()` |
| Brain-heart coupling (EEG-ECG) | `crossCorrelation()`, `waveletCoherence()` |
| Neurovascular coupling (EEG-fNIRS) | `coherence()`, `grangerCausality()` |
| EDA-EMG sympathetic co-activation | `slidingCrossCorrelation()`, `couplingMatrix()` |
| Movement-brain synchrony (MoCap-EEG) | `waveletPLV()`, `waveletCoherence()` |

## Dependencies

- **R** (>= 4.2)
- **PhysioCore** -- core data structures and accessors
- **methods**, **stats** -- base R infrastructure
- **SummarizedExperiment**, **S4Vectors** -- Bioconductor infrastructure

Optional (in Suggests):

- **signal** -- additional DSP functions
- **PhysioPreprocess** -- signal preprocessing
- **ggplot2** -- visualization
- **knitr**, **rmarkdown** -- vignettes

## PhysioExperiment Ecosystem

PhysioCrossModal is part of the PhysioExperiment ecosystem, a suite of R packages for multi-modal physiological signal analysis:

| Package | Description |
|---------|-------------|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| [PhysioPreprocess](https://github.com/x-biosignal/PhysioPreprocess) | Preprocessing (filters, ICA, resampling) |
| [PhysioAnalysis](https://github.com/x-biosignal/PhysioAnalysis) | Analysis and visualization |
| **PhysioCrossModal** | Cross-modal coupling and connectivity |
| [PhysioMoCap](https://github.com/x-biosignal/PhysioMoCap) | Motion capture and biomechanics |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to browse all available packages.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Author

Yusuke Matsui

## Governance & support

Part of the [Physio ecosystem](https://x-biosignal.r-universe.dev). Community and
policy documents live in the umbrella repository:

- [Code of Conduct](https://github.com/x-biosignal/PhysioExperiment/blob/main/CODE_OF_CONDUCT.md)
- [Contributing](https://github.com/x-biosignal/PhysioExperiment/blob/main/CONTRIBUTING.md)
- [Governance](https://github.com/x-biosignal/PhysioExperiment/blob/main/GOVERNANCE.md)
- [Support](https://github.com/x-biosignal/PhysioExperiment/blob/main/SUPPORT.md)
- [Security policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/SECURITY.md)
- [Deprecation & lifecycle policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/DEPRECATION.md)
