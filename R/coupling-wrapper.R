# Unified coupling analysis wrapper
#
# High-level dispatch function that routes to the appropriate coupling
# method based on a user-specified `method` argument.

#' Extract a scalar coupling statistic from a coupling result
#'
#' @param result List returned by a coupling function.
#' @param method Character method name used to produce the result.
#' @return Numeric scalar coupling statistic.
#' @noRd
.extract_coupling_stat <- function(result, method) {
  switch(method,
    coherence = max(abs(result$coherence)),
    multitaper_coherence = max(abs(result$coherence)),
    wavelet_coherence = mean(result$coherence),
    wavelet_plv = mean(result$plv),
    plv = result$plv,
    pli = result$pli,
    wpli = result$wpli,
    granger = result$gc_xy,
    crosscorrelation = abs(result$peak_correlation),
    dtf = result$dtf_xy,
    pdc = result$pdc_xy,
    transferentropy = result$te_xy,
    pac = result$pac,
    stop("Unknown method for stat extraction: '", method, "'", call. = FALSE)
  )
}

#' Unified interface for cross-modal coupling analysis
#'
#' Dispatches to the appropriate coupling function based on the `method`
#' parameter. Accepts the same flexible inputs as the underlying functions:
#' two numeric vectors (with `sr`), two PhysioExperiment objects, or a
#' MultiPhysioExperiment with named modalities.
#'
#' @param x Numeric vector, PhysioExperiment, or MultiPhysioExperiment.
#'   When `mpe` is provided this argument is ignored.
#' @param y Numeric vector or PhysioExperiment (NULL when `x` is an MPE or
#'   when `mpe` is provided).
#' @param mpe A \code{\link{MultiPhysioExperiment}} object. When supplied,
#'   `x` and `y` are ignored and signals are extracted from `mpe` using
#'   `modality_x` / `modality_y`.
#' @param modality_x,modality_y Character names of the modalities to extract
#'   from `mpe`.
#' @param channels_x,channels_y Integer channel indices to extract
#'   (default 1).
#' @param method Character string specifying the coupling method. One of:
#'   \code{"coherence"}, \code{"plv"}, \code{"pli"}, \code{"wpli"},
#'   \code{"granger"}, \code{"crosscorrelation"}, \code{"wavelet_coherence"},
#'   or \code{"wavelet_plv"}.
#' @param sr Numeric sampling rate in Hz. Required when `x` and `y` are
#'   numeric vectors.
#' @param ... Additional arguments passed to the specific coupling function
#'   (e.g. `freq_band`, `order`, `max_lag`, `nperseg`, etc.).
#' @param granger_method Named-only estimator selector used when
#'   `method = "granger"`: `"parametric"` (default) or `"nonparametric"`.
#'
#' @return The result from the dispatched coupling function. See individual
#'   function documentation for details:
#'   \itemize{
#'     \item \code{"coherence"}: see \code{\link{coherence}}
#'     \item \code{"plv"}: see \code{\link{phaseLockingValue}}
#'     \item \code{"pli"}: see \code{\link{phaseLagIndex}}
#'     \item \code{"wpli"}: see \code{\link{weightedPLI}}
#'     \item \code{"granger"}: see \code{\link{grangerCausality}}
#'     \item \code{"crosscorrelation"}: see \code{\link{crossCorrelation}}
#'     \item \code{"wavelet_coherence"}: see \code{\link{waveletCoherence}}
#'     \item \code{"wavelet_plv"}: see \code{\link{waveletPLV}}
#'   }
#'
#' @references
#' Carter, G. C. (1987). Coherence and time delay estimation.
#' \emph{Proceedings of the IEEE}, 75(2), 236--255.
#'
#' Lachaux, J.-P., Rodriguez, E., Martinerie, J., & Varela, F. J. (1999).
#' Measuring phase synchrony in brain signals. \emph{Human Brain Mapping},
#' 8(4), 194--208.
#'
#' Granger, C. W. J. (1969). Investigating causal relations by econometric
#' models and cross-spectral methods. \emph{Econometrica}, 37(3), 424--438.
#'
#' @seealso [coherence()], [phaseLockingValue()], [grangerCausality()],
#'   [crossCorrelation()], [surrogateTest()]
#' @export
#' @examples
#' # Numeric vectors
#' sr <- 500
#' t <- seq(0, 10, length.out = sr * 10)
#' x <- sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
#' y <- 0.8 * sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
#'
#' # Coherence
#' res <- couplingAnalysis(x, y, method = "coherence", sr = sr)
#'
#' # Cross-correlation
#' res <- couplingAnalysis(x, y, method = "crosscorrelation", sr = sr)
couplingAnalysis <- function(x, y = NULL, mpe = NULL,
                             modality_x = NULL, modality_y = NULL,
                             channels_x = 1L, channels_y = 1L,
                             method = c("coherence", "plv", "pli", "wpli",
                                        "granger", "crosscorrelation",
                                        "wavelet_coherence", "wavelet_plv",
                                        "multitaper_coherence", "dtf", "pdc",
                                        "transferentropy", "pac"),
                             sr = NULL, ...,
                             granger_method = c("parametric", "nonparametric")) {

  method <- match.arg(method)

  # Build the common arguments that every coupling function understands
  # If mpe is not explicitly provided but x is an MPE, treat it as the mpe path
  if (is.null(mpe) && inherits(x, "MultiPhysioExperiment")) {
    mpe <- x
    x <- NULL
  }

  if (!is.null(mpe)) {
    # MPE path: pass mpe as x, with modality names
    common_args <- list(
      x = mpe,
      y = NULL,
      modality_x = modality_x,
      modality_y = modality_y,
      channels_x = channels_x,
      channels_y = channels_y
    )
  } else {
    # x/y path: pass signals (numeric or PE) directly
    common_args <- list(
      x = x,
      y = y,
      channels_x = channels_x,
      channels_y = channels_y
    )
    # Only include sr if it was explicitly provided (numeric vectors need it;
    # PhysioExperiment objects carry their own)
    if (!is.null(sr)) {
      common_args$sr <- sr
    }
  }

  # Merge with user-supplied extra arguments
  extra_args <- list(...)

  # Dispatch to the appropriate function
  switch(method,
    coherence = {
      args <- c(common_args, extra_args)
      do.call(coherence, args)
    },
    plv = {
      args <- c(common_args, extra_args)
      do.call(phaseLockingValue, args)
    },
    pli = {
      args <- c(common_args, extra_args)
      do.call(phaseLagIndex, args)
    },
    wpli = {
      args <- c(common_args, extra_args)
      do.call(weightedPLI, args)
    },
    granger = {
      args <- c(common_args, extra_args)
      args$method <- match.arg(granger_method)
      do.call(grangerCausality, args)
    },
    crosscorrelation = {
      args <- c(common_args, extra_args)
      do.call(crossCorrelation, args)
    },
    wavelet_coherence = {
      args <- c(common_args, extra_args)
      do.call(waveletCoherence, args)
    },
    wavelet_plv = {
      args <- c(common_args, extra_args)
      do.call(waveletPLV, args)
    },
    multitaper_coherence = {
      args <- c(common_args, extra_args)
      do.call(multitaperCoherence, args)
    },
    dtf = {
      args <- c(common_args, extra_args)
      do.call(directedTransferFunction, args)
    },
    pdc = {
      args <- c(common_args, extra_args)
      do.call(partialDirectedCoherence, args)
    },
    transferentropy = {
      args <- c(common_args, extra_args)
      do.call(transferEntropy, args)
    },
    pac = {
      args <- c(common_args, extra_args)
      do.call(phaseAmplitudeCoupling, args)
    }
  )
}

#' Cached coupling analysis over a MultiPhysioExperiment
#'
#' Wraps \code{\link{couplingAnalysis}()} with a per-object memoisation cache
#' held in the \code{couplingResults} slot of a
#' \code{\link{MultiPhysioExperiment}}. Results are keyed by a
#' \code{\link[digest]{digest}} of the analysis arguments (method, modalities,
#' channels, and any extra arguments), so a repeated identical request returns
#' the stored result without recomputing.
#'
#' @param mpe A \code{MultiPhysioExperiment} object.
#' @param method Character coupling method (see \code{\link{couplingAnalysis}}).
#' @param modality_x,modality_y Character modality names to extract from `mpe`.
#' @param channels_x,channels_y Integer channel indices (default 1).
#' @param ... Additional arguments forwarded to \code{couplingAnalysis()} and
#'   folded into the cache key.
#' @param cache Logical; if \code{TRUE} (default) reuse and populate the cache.
#'   \code{FALSE} forces recomputation but still stores the fresh result.
#' @return A list with \code{result} (the coupling result), \code{mpe} (the
#'   object with its \code{couplingResults} cache updated) and \code{cached}
#'   (\code{TRUE} when the result was served from the cache).
#' @seealso [couplingAnalysis()], [couplingResults()]
#' @export
#' @examples
#' eeg <- PhysioExperiment(assays = list(raw = matrix(rnorm(1000), 500, 2)),
#'                         samplingRate = 250)
#' emg <- PhysioExperiment(assays = list(raw = matrix(rnorm(1000), 500, 2)),
#'                         samplingRate = 250)
#' mpe <- MultiPhysioExperiment(EEG = eeg, EMG = emg)
#' out  <- couplingAnalysisCached(mpe, "coherence", "EEG", "EMG")
#' out2 <- couplingAnalysisCached(out$mpe, "coherence", "EEG", "EMG")
#' out2$cached  # TRUE
couplingAnalysisCached <- function(mpe, method,
                                   modality_x = NULL, modality_y = NULL,
                                   channels_x = 1L, channels_y = 1L, ...,
                                   cache = TRUE) {
  stopifnot(inherits(mpe, "MultiPhysioExperiment"))
  key <- digest::digest(list(
    method = method, modality_x = modality_x, modality_y = modality_y,
    channels_x = channels_x, channels_y = channels_y, extra = list(...)
  ))
  store <- mpe@couplingResults
  if (cache && !is.null(store[[key]])) {
    return(list(result = store[[key]], mpe = mpe, cached = TRUE))
  }
  result <- couplingAnalysis(
    mpe = mpe, method = method,
    modality_x = modality_x, modality_y = modality_y,
    channels_x = channels_x, channels_y = channels_y, ...
  )
  store[[key]] <- result
  mpe@couplingResults <- store
  list(result = result, mpe = mpe, cached = FALSE)
}
