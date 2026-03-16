# Visualization functions for cross-modal coupling results
#
# All functions require ggplot2 (listed in Suggests) and return ggplot objects.

# Suppress R CMD check NOTEs for ggplot2 aesthetics
# These variables are column names used in aes() calls and are not global
# bindings.
utils::globalVariables(c("Col", "Row", "Value", "Frequency", "Coherence",
                          "Time", "Correlation",
                          "WavCoh_Value", "WavCoh_Freq", "WavCoh_Time",
                          "COI_Freq"))

#' Plot a coupling matrix as a heatmap
#'
#' Displays a matrix of coupling values (e.g., coherence between all channel
#' pairs) as a colour-coded heatmap using \code{ggplot2::geom_tile()}.
#'
#' @param mat Numeric matrix of coupling values. Row and column names, if
#'   present, are used as axis labels.
#' @param title Character string for the plot title
#'   (default \code{"Coupling Matrix"}).
#' @param low_colour Character colour for the low end of the scale
#'   (default \code{"white"}).
#' @param high_colour Character colour for the high end of the scale
#'   (default \code{"#2166AC"}).
#' @param ... Additional arguments passed to \code{ggplot2::theme()}.
#'
#' @references
#' Wickham, H. (2016). \emph{ggplot2: Elegant Graphics for Data Analysis}.
#' Springer-Verlag New York. \doi{10.1007/978-3-319-24277-4}
#'
#' @seealso [coherenceMatrix()], [couplingMatrix()],
#'   [surrogateMatrixTest()]
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' mat <- matrix(runif(16), nrow = 4,
#'               dimnames = list(paste0("Ch", 1:4), paste0("Ch", 1:4)))
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotCouplingMatrix(mat)
#' }
plotCouplingMatrix <- function(mat, title = "Coupling Matrix",
                               low_colour = "white",
                               high_colour = "#2166AC", ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("'mat' must be a numeric matrix", call. = FALSE)
  }

  # Use row/col names or generate defaults
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("R", seq_len(nrow(mat)))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("C", seq_len(ncol(mat)))
  }

  # Convert to long format for ggplot
  df <- expand.grid(
    Row = factor(rownames(mat), levels = rev(rownames(mat))),
    Col = factor(colnames(mat), levels = colnames(mat))
  )
  df$Value <- as.vector(t(mat[rev(seq_len(nrow(mat))), , drop = FALSE]))

  ggplot2::ggplot(df, ggplot2::aes(x = Col, y = Row, fill = Value)) +
    ggplot2::geom_tile(colour = "grey80") +
    ggplot2::scale_fill_gradient(low = low_colour, high = high_colour) +
    ggplot2::labs(title = title, x = "", y = "", fill = "Coupling") +
    ggplot2::theme_minimal(...) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}


#' Plot coherence spectrum
#'
#' Displays the coherence as a function of frequency from the result of
#' \code{\link{coherence}}. Optionally draws the 95\% confidence limit as
#' a horizontal dashed line.
#'
#' @param result A list returned by \code{\link{coherence}}, containing at
#'   least \code{$coherence} and \code{$frequencies}. If
#'   \code{$confidence_limit} is present it is drawn as a threshold line.
#' @param show_threshold Logical; if TRUE (default) and a confidence limit
#'   is available, draw it on the plot.
#' @param title Character string for the plot title
#'   (default \code{"Coherence Spectrum"}).
#' @param ... Additional arguments (currently unused).
#'
#' @references
#' Wickham, H. (2016). \emph{ggplot2: Elegant Graphics for Data Analysis}.
#' Springer-Verlag New York. \doi{10.1007/978-3-319-24277-4}
#'
#' @seealso [coherence()], [multitaperCoherence()],
#'   [plotCouplingMatrix()]
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' sr <- 500
#' t <- seq(0, 10, length.out = sr * 10)
#' x <- sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
#' y <- 0.8 * sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
#' res <- coherence(x, y, sr = sr)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotCoherenceSpectrum(res)
#' }
plotCoherenceSpectrum <- function(result, show_threshold = TRUE,
                                  title = "Coherence Spectrum", ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  if (!is.list(result) || is.null(result$coherence) ||
      is.null(result$frequencies)) {
    stop("'result' must be a list with 'coherence' and 'frequencies' components",
         call. = FALSE)
  }

  df <- data.frame(
    Frequency = result$frequencies,
    Coherence = result$coherence
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Frequency, y = Coherence)) +
    ggplot2::geom_line(colour = "#2166AC", linewidth = 0.7) +
    ggplot2::labs(title = title, x = "Frequency (Hz)",
                  y = "Coherence") +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal()

  # Add confidence threshold line
  if (show_threshold && !is.null(result$confidence_limit)) {
    p <- p + ggplot2::geom_hline(
      yintercept = result$confidence_limit,
      linetype = "dashed", colour = "red", linewidth = 0.5
    )
  }

  p
}


#' Plot time-varying coupling from sliding-window analysis
#'
#' Displays the peak cross-correlation over time from the result of
#' \code{\link{slidingCrossCorrelation}}.
#'
#' @param result A list returned by \code{\link{slidingCrossCorrelation}},
#'   containing at least \code{$times} and \code{$peak_correlations}.
#' @param title Character string for the plot title
#'   (default \code{"Coupling Time Course"}).
#' @param ... Additional arguments (currently unused).
#'
#' @references
#' Wickham, H. (2016). \emph{ggplot2: Elegant Graphics for Data Analysis}.
#' Springer-Verlag New York. \doi{10.1007/978-3-319-24277-4}
#'
#' @seealso [slidingCrossCorrelation()], [crossCorrelation()],
#'   [plotCoherenceSpectrum()]
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' sr <- 500
#' set.seed(1)
#' x <- rnorm(5000)
#' y <- c(rep(0, 10), x[1:4990])
#' res <- slidingCrossCorrelation(x, y, sr = sr,
#'                                 window_sec = 1, step_sec = 0.5)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotCouplingTimecourse(res)
#' }
plotCouplingTimecourse <- function(result,
                                   title = "Coupling Time Course", ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  if (!is.list(result) || is.null(result$times) ||
      is.null(result$peak_correlations)) {
    stop("'result' must be a list with 'times' and 'peak_correlations' components",
         call. = FALSE)
  }

  df <- data.frame(
    Time = result$times,
    Correlation = result$peak_correlations
  )

  ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Correlation)) +
    ggplot2::geom_line(colour = "#2166AC", linewidth = 0.7) +
    ggplot2::labs(title = title, x = "Time (s)",
                  y = "Peak Correlation") +
    ggplot2::theme_minimal()
}


#' Plot wavelet coherence time-frequency map
#'
#' Displays wavelet coherence (or wavelet PLV) as a filled time-frequency
#' heatmap with optional Cone of Influence (COI) overlay.
#'
#' @param result List returned by \code{\link{waveletCoherence}} or
#'   \code{\link{waveletPLV}}.
#' @param show_coi Logical; overlay COI boundary (default \code{TRUE} if
#'   \code{coi} is present in the result).
#' @param title Character plot title (default \code{"Wavelet Coherence"}).
#' @param fill_label Character legend label (default \code{"Coherence"}).
#' @param ... Additional arguments passed to \code{ggplot2::theme()}.
#'
#' @references
#' Wickham, H. (2016). \emph{ggplot2: Elegant Graphics for Data Analysis}.
#' Springer-Verlag New York. \doi{10.1007/978-3-319-24277-4}
#'
#' @seealso [waveletCoherence()], [waveletPLV()],
#'   [plotCoherenceSpectrum()]
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' sr <- 200
#' t <- seq(0, 2, length.out = sr * 2)
#' x <- sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
#' y <- 0.8 * sin(2 * pi * 10 * t) + 0.3 * rnorm(length(t))
#' result <- waveletCoherence(x, y, sr = sr, frequencies = seq(5, 20))
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotWaveletCoherence(result)
#' }
plotWaveletCoherence <- function(result,
                                  show_coi = TRUE,
                                  title = "Wavelet Coherence",
                                  fill_label = "Coherence",
                                  ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  # Auto-detect PLV vs coherence
  if (!is.null(result$plv) && is.null(result$coherence)) {
    mat <- result$plv
    if (fill_label == "Coherence") fill_label <- "PLV"
  } else if (!is.null(result$coherence)) {
    mat <- result$coherence
  } else {
    stop("'result' must contain a 'coherence' or 'plv' matrix", call. = FALSE)
  }

  if (is.null(result$frequencies) || is.null(result$times)) {
    stop("'result' must contain 'frequencies' and 'times' vectors",
         call. = FALSE)
  }

  # Convert to long-format data.frame
  n_time <- nrow(mat)
  n_freq <- ncol(mat)
  df <- data.frame(
    WavCoh_Time = rep(result$times, times = n_freq),
    WavCoh_Freq = rep(result$frequencies, each = n_time),
    WavCoh_Value = as.vector(mat)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = WavCoh_Time, y = WavCoh_Freq,
                                         fill = WavCoh_Value)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_viridis_c(option = "inferno", limits = c(0, 1)) +
    ggplot2::labs(title = title, x = "Time (s)", y = "Frequency (Hz)",
                  fill = fill_label) +
    ggplot2::theme_minimal(...)

  # Add COI overlay
  if (show_coi && !is.null(result$coi)) {
    coi_df <- data.frame(
      WavCoh_Time = result$times,
      COI_Freq = result$coi
    )
    # Clamp COI to frequency range for display
    freq_max <- max(result$frequencies)
    coi_df$COI_Freq <- pmin(coi_df$COI_Freq, freq_max)

    p <- p + ggplot2::geom_line(
      data = coi_df,
      ggplot2::aes(x = WavCoh_Time, y = COI_Freq),
      inherit.aes = FALSE,
      colour = "white", linewidth = 0.8, linetype = "dashed"
    )
  }

  p
}
