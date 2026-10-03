#' PhysioCrossModal: Cross-Modal Coupling Analysis for Physiological Signals
#'
#' Coupling, connectivity, synchrony and multi-block fusion between physiological
#' signals of different modalities (EEG, EMG, ECG, EDA, MoCap, fNIRS, ...). It
#' adds the `MultiPhysioExperiment` container for holding several
#' `PhysioExperiment` objects recorded at different sampling rates with temporal
#' alignment, and an analysis surface for spectral, phase, directed and
#' time-domain coupling plus surrogate-based inference.
#'
#' @section Container and alignment:
#' Build and inspect with [MultiPhysioExperiment()], [experiments()],
#' [modalities()], [nModalities()], [samplingRates()], [alignment()],
#' [mergePhysio()]. Put modalities on a common clock with [alignSignals()],
#' [alignToRate()], [elasticAlign()], [warpApply()].
#'
#' @section Spectral coupling:
#' [coherence()], [coherenceMatrix()], [crossSpectrum()],
#' [multitaperCoherence()], [waveletCoherence()], [laggedCoherence()],
#' [crossWaveletTransform()].
#'
#' @section Phase synchrony:
#' [phaseLockingValue()], [phaseLagIndex()], [weightedPLI()], [ciPLV()],
#' [waveletPLV()], [pairwisePhaseConsistency()], [phaseSlopeIndex()].
#'
#' @section Phase-amplitude coupling:
#' [phaseAmplitudeCoupling()], [comodulogram()], [modulationIndex()].
#'
#' @section Directed and time-domain coupling:
#' [grangerCausality()], [partialDirectedCoherence()],
#' [directedTransferFunction()], [transferEntropy()], [orthogonalizedAEC()];
#' [crossCorrelation()], [slidingCrossCorrelation()], [mutualInformation()].
#'
#' @section Unified coupling API:
#' [couplingAnalysis()] dispatches to any single method by name;
#' [couplingMatrix()] computes a pairwise matrix for any supported method.
#'
#' @section Multi-block fusion:
#' [jive()], [multipleFactorAnalysis()], [cca()], [plsBlocks()], [coupledNMF()],
#' [fuseBlocks()], [crossBlockFactor()], [rvCoefficient()],
#' [distanceCorrelation()], [representationalSimilarity()],
#' [multimodalGaitIndex()].
#'
#' @section Statistical inference and visualization:
#' [surrogateTest()], [surrogateMatrixTest()], [bootstrapCI()],
#' [lodoGeneralization()]; [plotCoherenceSpectrum()], [plotCouplingMatrix()],
#' [plotCouplingTimecourse()], [plotWaveletCoherence()], [plotComodulogram()].
#'
#' @section Simulated data for trying the methods:
#' [make_coupled_signals()], [make_directed_signals()], [make_eeg_emg()].
#'
#' @section Where to go next:
#' The data model and single-modality analysis live in \pkg{PhysioExperiment}
#' and \pkg{PhysioAnalysis}. Coupling measures are associations, not proof of a
#' biological mechanism: read each function's help for its interpretation
#' caveats (e.g. PAC and Granger causality). See
#' `vignette("introduction", package = "PhysioCrossModal")` and
#' `vignette("coupling-analysis", package = "PhysioCrossModal")`.
#'
#' @keywords internal
"_PACKAGE"
