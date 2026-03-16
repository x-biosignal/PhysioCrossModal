# Multi-channel coupling matrices
#
# Functions for computing coupling between all pairs of channels,
# returning matrix representations suitable for plotting with
# plotCouplingMatrix().

#' Coherence matrix across channel pairs
#'
#' Computes magnitude-squared coherence for all pairs of channels between
#' two sets of signals (or between all channels within a single set).
#' Returns a matrix of peak coherence values plus the full spectra.
#'
#' @param x A PhysioExperiment or MultiPhysioExperiment object.
#' @param y A PhysioExperiment, or NULL when \code{x} is an MPE.
#' @param sr Numeric sampling rate in Hz (used only for numeric inputs).
#' @param channels_x Integer vector of channel indices for x, or NULL
#'   for all channels.
#' @param channels_y Integer vector of channel indices for y, or NULL
#'   for all channels.
#' @param modality_x,modality_y Character modality names for MPE input.
#' @param ... Additional arguments passed to \code{\link{coherence}}
#'   (e.g. \code{nperseg}, \code{freq_range}).
#'
#' @return A list with components:
#'   \describe{
#'     \item{matrix}{Numeric matrix of peak coherence values, with
#'       dimensions \code{[n_channels_x x n_channels_y]}.}
#'     \item{spectra}{List of lists containing the full coherence result
#'       for each pair.}
#'     \item{frequencies}{Numeric vector of frequencies from the first pair.}
#'     \item{channel_names_x}{Character vector of x channel names.}
#'     \item{channel_names_y}{Character vector of y channel names.}
#'   }
#'
#' @references
#' Carter, G. C. (1987). Coherence and time delay estimation.
#' \emph{Proceedings of the IEEE}, 75(2), 236--255.
#'
#' @seealso [coherence()], [couplingMatrix()], [plotCouplingMatrix()],
#'   [surrogateMatrixTest()]
#' @export
#' @examples
#' sr <- 200; n <- sr * 2
#' pe1 <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
#'   samplingRate = sr
#' )
#' pe2 <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
#'   samplingRate = sr
#' )
#' result <- coherenceMatrix(pe1, pe2, nperseg = 64L)
#' result$matrix
coherenceMatrix <- function(x, y = NULL, sr = NULL,
                            channels_x = NULL, channels_y = NULL,
                            modality_x = NULL, modality_y = NULL,
                            ...) {

  # Resolve input to two PE objects
  pe_pair <- .resolve_pe_pair(x, y, modality_x, modality_y)
  pe_x <- pe_pair$pe_x
  pe_y <- pe_pair$pe_y

  # Resolve channel indices
  n_ch_x <- ncol(SummarizedExperiment::assay(pe_x, 1L))
  n_ch_y <- ncol(SummarizedExperiment::assay(pe_y, 1L))

  if (is.null(channels_x)) channels_x <- seq_len(n_ch_x)
  if (is.null(channels_y)) channels_y <- seq_len(n_ch_y)

  # Channel names
  cd_x <- SummarizedExperiment::colData(pe_x)
  cd_y <- SummarizedExperiment::colData(pe_y)
  ch_names_x <- if ("label" %in% colnames(cd_x)) {
    as.character(cd_x$label[channels_x])
  } else {
    paste0("x_ch", channels_x)
  }
  ch_names_y <- if ("label" %in% colnames(cd_y)) {
    as.character(cd_y$label[channels_y])
  } else {
    paste0("y_ch", channels_y)
  }

  # Compute coherence for all pairs
  mat <- matrix(NA_real_, nrow = length(channels_x), ncol = length(channels_y),
                dimnames = list(ch_names_x, ch_names_y))
  spectra <- vector("list", length(channels_x))
  names(spectra) <- ch_names_x
  freqs <- NULL

  for (i in seq_along(channels_x)) {
    spectra[[i]] <- vector("list", length(channels_y))
    names(spectra[[i]]) <- ch_names_y
    for (j in seq_along(channels_y)) {
      result <- coherence(pe_x, pe_y,
                          channels_x = channels_x[i],
                          channels_y = channels_y[j], ...)
      mat[i, j] <- max(result$coherence)
      spectra[[i]][[j]] <- result
      if (is.null(freqs)) freqs <- result$frequencies
    }
  }

  list(
    matrix = mat,
    spectra = spectra,
    frequencies = freqs,
    channel_names_x = ch_names_x,
    channel_names_y = ch_names_y
  )
}


#' Generic coupling matrix across channel pairs
#'
#' Computes any coupling method for all pairs of channels between two
#' sets of signals, extracting a scalar statistic for each pair. This
#' is the multi-channel generalisation of \code{\link{couplingAnalysis}}.
#'
#' @param x A PhysioExperiment or MultiPhysioExperiment object.
#' @param y A PhysioExperiment, or NULL when \code{x} is an MPE.
#' @param sr Numeric sampling rate in Hz.
#' @param method Character coupling method, one of \code{"coherence"},
#'   \code{"plv"}, \code{"pli"}, \code{"wpli"}, \code{"granger"},
#'   \code{"crosscorrelation"}.
#' @param channels_x,channels_y Integer vectors of channel indices, or
#'   NULL for all channels.
#' @param modality_x,modality_y Character modality names for MPE input.
#' @param ... Additional arguments passed to the coupling function.
#'
#' @return A list with components:
#'   \describe{
#'     \item{matrix}{Numeric matrix of coupling values.}
#'     \item{method}{Character method used.}
#'     \item{channel_names_x}{Character vector of x channel names.}
#'     \item{channel_names_y}{Character vector of y channel names.}
#'   }
#'
#' @references
#' Carter, G. C. (1987). Coherence and time delay estimation.
#' \emph{Proceedings of the IEEE}, 75(2), 236--255.
#'
#' @seealso [couplingAnalysis()], [coherenceMatrix()],
#'   [plotCouplingMatrix()], [surrogateMatrixTest()]
#' @export
#' @examples
#' sr <- 200; n <- sr * 2
#' pe1 <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
#'   samplingRate = sr
#' )
#' pe2 <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(n * 2), nrow = n)),
#'   samplingRate = sr
#' )
#' result <- couplingMatrix(pe1, pe2, method = "coherence", nperseg = 64L)
#' result$matrix
couplingMatrix <- function(x, y = NULL, sr = NULL, method,
                           channels_x = NULL, channels_y = NULL,
                           modality_x = NULL, modality_y = NULL,
                           ...) {

  # Resolve input to two PE objects
  pe_pair <- .resolve_pe_pair(x, y, modality_x, modality_y)
  pe_x <- pe_pair$pe_x
  pe_y <- pe_pair$pe_y

  # Resolve channel indices
  n_ch_x <- ncol(SummarizedExperiment::assay(pe_x, 1L))
  n_ch_y <- ncol(SummarizedExperiment::assay(pe_y, 1L))

  if (is.null(channels_x)) channels_x <- seq_len(n_ch_x)
  if (is.null(channels_y)) channels_y <- seq_len(n_ch_y)

  # Channel names
  cd_x <- SummarizedExperiment::colData(pe_x)
  cd_y <- SummarizedExperiment::colData(pe_y)
  ch_names_x <- if ("label" %in% colnames(cd_x)) {
    as.character(cd_x$label[channels_x])
  } else {
    paste0("x_ch", channels_x)
  }
  ch_names_y <- if ("label" %in% colnames(cd_y)) {
    as.character(cd_y$label[channels_y])
  } else {
    paste0("y_ch", channels_y)
  }

  mat <- matrix(NA_real_, nrow = length(channels_x), ncol = length(channels_y),
                dimnames = list(ch_names_x, ch_names_y))

  for (i in seq_along(channels_x)) {
    for (j in seq_along(channels_y)) {
      result <- couplingAnalysis(pe_x, pe_y,
                                 channels_x = channels_x[i],
                                 channels_y = channels_y[j],
                                 method = method, ...)
      mat[i, j] <- .extract_coupling_stat(result, method)
    }
  }

  list(
    matrix = mat,
    method = method,
    channel_names_x = ch_names_x,
    channel_names_y = ch_names_y
  )
}


# ---- Internal helpers --------------------------------------------------------

#' Resolve input to a pair of PhysioExperiment objects
#'
#' @param x PhysioExperiment or MultiPhysioExperiment.
#' @param y PhysioExperiment or NULL.
#' @param modality_x Character modality name.
#' @param modality_y Character modality name.
#' @return Named list with \code{pe_x} and \code{pe_y}.
#' @noRd
.resolve_pe_pair <- function(x, y, modality_x, modality_y) {
  if (inherits(x, "MultiPhysioExperiment")) {
    if (is.null(modality_x) || is.null(modality_y)) {
      stop("modality_x and modality_y must be specified for MultiPhysioExperiment",
           call. = FALSE)
    }
    exps <- experiments(x)
    if (!modality_x %in% names(exps)) {
      stop("modality_x '", modality_x, "' not found", call. = FALSE)
    }
    if (!modality_y %in% names(exps)) {
      stop("modality_y '", modality_y, "' not found", call. = FALSE)
    }
    return(list(pe_x = exps[[modality_x]], pe_y = exps[[modality_y]]))
  }

  if (inherits(x, "PhysioExperiment")) {
    if (is.null(y) || !inherits(y, "PhysioExperiment")) {
      stop("When x is a PhysioExperiment, y must also be a PhysioExperiment",
           call. = FALSE)
    }
    return(list(pe_x = x, pe_y = y))
  }

  stop("x must be a PhysioExperiment or MultiPhysioExperiment", call. = FALSE)
}
