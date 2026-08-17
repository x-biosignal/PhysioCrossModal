#' MultiPhysioExperiment class definition
#'
#' The `MultiPhysioExperiment` class holds multiple `PhysioExperiment` objects
#' representing different signal modalities (e.g., EEG, EMG, ECG) recorded
#' simultaneously but potentially at different sampling rates. It provides
#' temporal alignment metadata and a cache for coupling analysis results.
#'
#' @slot experiments Named list of \code{PhysioExperiment} objects.
#' @slot alignment \code{\link[S4Vectors]{DataFrame}} with temporal alignment
#'   metadata (one row per modality).
#' @slot sampleMap \code{\link[S4Vectors]{DataFrame}} mapping samples across
#'   modalities.
#' @slot couplingResults List of cached coupling analysis results.
#'
#' @references
#' Huber, W., et al. (2015). "Orchestrating high-throughput genomic analysis
#' with Bioconductor." Nature Methods, 12(2), 115-121.
#' doi:10.1038/nmeth.3252
#'
#' @seealso [MultiPhysioExperiment()], [experiments()], [modalities()],
#'   [alignSignals()]
#' @exportClass MultiPhysioExperiment
setClass(
  "MultiPhysioExperiment",
  slots = c(
    experiments = "list",
    alignment = "DataFrame",
    sampleMap = "DataFrame",
    couplingResults = "list"
  ),
  prototype = list(
    experiments = list(),
    alignment = S4Vectors::DataFrame(),
    sampleMap = S4Vectors::DataFrame(),
    couplingResults = list()
  ),
  validity = function(object) {
    msgs <- character()
    exps <- object@experiments

    if (length(exps) == 0) {
      msgs <- c(msgs, "experiments must contain at least one PhysioExperiment")
    }

    if (!all(vapply(exps, inherits, logical(1), "PhysioExperiment"))) {
      msgs <- c(msgs, "all experiments must be PhysioExperiment objects")
    }

    if (is.null(names(exps)) || any(names(exps) == "")) {
      msgs <- c(msgs, "all experiments must be named")
    }

    if (length(msgs) == 0) TRUE else msgs
  }
)

# ---- Constructor ---------------------------------------------------------

#' Construct a MultiPhysioExperiment object
#'
#' Creates a container that holds multiple \code{PhysioExperiment} objects
#' recorded simultaneously at potentially different sampling rates.
#'
#' @param ... Named \code{PhysioExperiment} objects, one per modality.
#' @param experiments Alternatively, a named list of \code{PhysioExperiment}
#'   objects.  If both \code{...} and \code{experiments} are provided, they
#'   are combined.
#' @param alignment Optional \code{\link[S4Vectors]{DataFrame}} with temporal
#'   alignment metadata.  When \code{NULL} (the default), a default alignment
#'   table is built from the supplied experiments.
#' @return A \code{\link{MultiPhysioExperiment-class}} instance containing
#'   the supplied experiments, alignment metadata, and an empty coupling
#'   results cache.
#'
#' @references
#' Huber, W., et al. (2015). "Orchestrating high-throughput genomic analysis
#' with Bioconductor." Nature Methods, 12(2), 115-121.
#' doi:10.1038/nmeth.3252
#'
#' @seealso [experiments()], [modalities()], [samplingRates()],
#'   [alignSignals()], [couplingAnalysis()]
#' @export
#' @examples
#' # Create two PhysioExperiment objects with different sampling rates
#' eeg_data <- matrix(rnorm(500 * 4), nrow = 500, ncol = 4)
#' emg_data <- matrix(rnorm(1000 * 2), nrow = 1000, ncol = 2)
#'
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   colData = S4Vectors::DataFrame(
#'     label = c("Fz", "Cz", "Pz", "Oz"),
#'     type = rep("EEG", 4)
#'   ),
#'   samplingRate = 250
#' )
#'
#' pe_emg <- PhysioExperiment(
#'   assays = list(raw = emg_data),
#'   colData = S4Vectors::DataFrame(
#'     label = c("EMG1", "EMG2"),
#'     type = rep("EMG", 2)
#'   ),
#'   samplingRate = 1000
#' )
#'
#' # Construct MultiPhysioExperiment
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)
#' mpe
MultiPhysioExperiment <- function(..., experiments = list(),
                                   alignment = NULL) {
  dots <- list(...)

  # Combine dots and explicit experiments list
  all_exps <- c(dots, experiments)

  # Build default alignment from experiments
  if (is.null(alignment)) {
    alignment <- .build_default_alignment(all_exps)
  }

  methods::new(
    "MultiPhysioExperiment",
    experiments = all_exps,
    alignment = alignment,
    sampleMap = S4Vectors::DataFrame(),
    couplingResults = list()
  )
}

#' Build default alignment DataFrame from experiments
#'
#' @param exps Named list of PhysioExperiment objects.
#' @return A DataFrame with one row per modality.
#' @keywords internal
.build_default_alignment <- function(exps) {
  if (length(exps) == 0) {
    return(S4Vectors::DataFrame())
  }
  nms <- names(exps)
  if (is.null(nms)) nms <- rep("", length(exps))
  srs <- vapply(exps, function(e) {
    if (inherits(e, "PhysioExperiment")) samplingRate(e) else NA_real_
  }, numeric(1))
  S4Vectors::DataFrame(
    modality = nms,
    samplingRate = srs,
    offset = rep(0, length(exps))
  )
}

# ---- Generics and Accessors ---------------------------------------------

#' Get or set the experiments list
#'
#' @param x A \code{MultiPhysioExperiment} object.
#' @param value A named list of \code{PhysioExperiment} objects.
#' @return For the getter, a named list of \code{PhysioExperiment} objects.
#'   For the setter, the modified \code{MultiPhysioExperiment} (returned
#'   invisibly).
#'
#' @seealso [MultiPhysioExperiment()], [modalities()], [samplingRates()]
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   samplingRate = 250
#' )
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg)
#' experiments(mpe)
setGeneric("experiments", function(x) standardGeneric("experiments"))

#' @rdname experiments
#' @export
setMethod("experiments", "MultiPhysioExperiment", function(x) x@experiments)

#' @rdname experiments
#' @export
setGeneric("experiments<-", function(x, value) standardGeneric("experiments<-"))

#' @rdname experiments
#' @export
setReplaceMethod("experiments", "MultiPhysioExperiment", function(x, value) {
  x@experiments <- value
  methods::validObject(x)
  x
})

#' Get modality names
#'
#' Returns a character vector of the names of the modalities stored in
#' the object.
#'
#' @param x A \code{MultiPhysioExperiment} object.
#' @return Character vector of modality names (e.g., \code{c("EEG", "EMG")}).
#'
#' @seealso [experiments()], [nModalities()], [samplingRates()]
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
#' emg_data <- matrix(rnorm(1000 * 2), nrow = 1000, ncol = 2)
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   samplingRate = 250
#' )
#' pe_emg <- PhysioExperiment(
#'   assays = list(raw = emg_data),
#'   samplingRate = 1000
#' )
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)
#' modalities(mpe)
setGeneric("modalities", function(x) standardGeneric("modalities"))

#' @rdname modalities
#' @export
setMethod("modalities", "MultiPhysioExperiment", function(x) {
  names(x@experiments)
})

#' Get sampling rates for all modalities
#'
#' Returns a named numeric vector of sampling rates, one per modality.
#'
#' @param x A \code{MultiPhysioExperiment} object.
#' @return Named numeric vector of sampling rates in Hz (e.g.,
#'   \code{c(EEG = 500, EMG = 1000)}).
#'
#' @seealso [modalities()], [experiments()], [alignToRate()]
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
#' emg_data <- matrix(rnorm(1000 * 2), nrow = 1000, ncol = 2)
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   samplingRate = 250
#' )
#' pe_emg <- PhysioExperiment(
#'   assays = list(raw = emg_data),
#'   samplingRate = 1000
#' )
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)
#' samplingRates(mpe)
setGeneric("samplingRates", function(x) standardGeneric("samplingRates"))

#' @rdname samplingRates
#' @export
setMethod("samplingRates", "MultiPhysioExperiment", function(x) {
  exps <- x@experiments
  srs <- vapply(exps, samplingRate, numeric(1))
  names(srs) <- names(exps)
  srs
})

#' Get number of modalities
#'
#' @param x A \code{MultiPhysioExperiment} object.
#' @return Integer scalar: the number of modalities stored in the object.
#'
#' @seealso [modalities()], [experiments()]
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   samplingRate = 250
#' )
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg)
#' nModalities(mpe)
setGeneric("nModalities", function(x) standardGeneric("nModalities"))

#' @rdname nModalities
#' @export
setMethod("nModalities", "MultiPhysioExperiment", function(x) {
  length(x@experiments)
})

#' Get or set temporal alignment metadata
#'
#' @param x A \code{MultiPhysioExperiment} object.
#' @param value A \code{\link[S4Vectors]{DataFrame}} with alignment metadata.
#' @return A \code{\link[S4Vectors]{DataFrame}}.
#' @seealso [alignSignals()], [experiments()], [samplingRates()]
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   samplingRate = 250
#' )
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg)
#' alignment(mpe)
setGeneric("alignment", function(x) standardGeneric("alignment"))

#' @rdname alignment
#' @export
setMethod("alignment", "MultiPhysioExperiment", function(x) x@alignment)

#' @rdname alignment
#' @export
setGeneric("alignment<-", function(x, value) standardGeneric("alignment<-"))

#' @rdname alignment
#' @export
setReplaceMethod("alignment", "MultiPhysioExperiment", function(x, value) {
  x@alignment <- value
  x
})

# ---- Show method ---------------------------------------------------------

#' Show method for MultiPhysioExperiment
#'
#' Displays a human-readable summary of the object.
#'
#' @param object A \code{MultiPhysioExperiment} object.
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   samplingRate = 250
#' )
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg)
#' mpe
setMethod("show", "MultiPhysioExperiment", function(object) {
  cat("class: MultiPhysioExperiment\n")
  n <- nModalities(object)
  cat("modalities(", n, "): ",
      paste(utils::head(modalities(object), 5), collapse = ", "),
      if (n > 5) " ..." else "", "\n", sep = "")

  srs <- samplingRates(object)
  cat("samplingRates: ",
      paste(paste0(names(srs), "=", srs, "Hz"), collapse = ", "),
      "\n", sep = "")

  exps <- experiments(object)
  for (nm in names(exps)) {
    pe <- exps[[nm]]
    n_ch <- ncol(SummarizedExperiment::assay(pe))
    n_tp <- nrow(SummarizedExperiment::assay(pe))
    cat("  ", nm, ": ", n_tp, " timepoints x ", n_ch, " channels\n", sep = "")
  }

  # Coupling results summary
  cr <- object@couplingResults
  if (length(cr) > 0) {
    cat("couplingResults(", length(cr), "): ",
        paste(utils::head(names(cr), 3), collapse = ", "),
        if (length(cr) > 3) " ..." else "", "\n", sep = "")
  }
})

# ---- Subsetting ----------------------------------------------------------

#' Extract a single modality from a MultiPhysioExperiment
#'
#' @param x A \code{MultiPhysioExperiment} object.
#' @param i Modality name (character) or index (integer).
#' @param j Not used.
#' @param ... Additional arguments (not used).
#' @return A \code{PhysioExperiment} object.
#' @export
#' @examples
#' eeg_data <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
#' pe_eeg <- PhysioExperiment(
#'   assays = list(raw = eeg_data),
#'   samplingRate = 250
#' )
#' mpe <- MultiPhysioExperiment(EEG = pe_eeg)
#' mpe[["EEG"]]
setMethod("[[", "MultiPhysioExperiment", function(x, i, j, ...) {
  x@experiments[[i]]
})
