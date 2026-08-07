# Elastic (SRVF) curve registration: separate phase (warping) from amplitude by
# working in the Square-Root Velocity Function representation (Srivastava et al.
# 2011). A pure-R dynamic-programming warping is provided; when the fdasrvf
# package is installed it is used for higher-accuracy alignment.

# Fast piecewise-linear interpolation of (t, y) at xout.
.lin_interp <- function(t, y, xout) {
  idx <- findInterval(xout, t, all.inside = TRUE)
  w <- (xout - t[idx]) / (t[idx + 1] - t[idx])
  y[idx] * (1 - w) + y[idx + 1] * w
}

# Square-Root Velocity Function of a curve f sampled on grid t.
.srvf <- function(f, t) {
  d <- diff(f) / diff(t)
  d <- c(d, d[length(d)])
  sign(d) * sqrt(abs(d))
}

# Reconstruct a curve from its SRVF: f'(t) = q(t)|q(t)|.
.srvf_to_curve <- function(q, t, f0 = 0) {
  fp <- q * abs(q)
  f0 + c(0, cumsum((fp[-1] + fp[-length(fp)]) / 2 * diff(t)))
}

# Action of a warping gamma on a SRVF: (q o gamma) * sqrt(gamma').
.warp_srvf <- function(q, gamma, t) {
  gp <- diff(gamma) / diff(t); gp <- c(gp, gp[length(gp)]); gp[gp < 0] <- 0
  .lin_interp(t, q, gamma) * sqrt(gp)
}

# Dynamic-programming optimal warping aligning q2 to q1 on grid t. Returns a
# monotone, boundary-preserving gamma on the grid.
.dp_warp <- function(q1, q2, t, max_step = 6L) {
  N <- length(t)
  E <- matrix(Inf, N, N); E[1, 1] <- 0
  Pi <- matrix(0L, N, N); Pj <- matrix(0L, N, N)
  for (i in 2:N) for (j in 2:N) {
    best <- Inf; bk <- 0L; bl <- 0L
    for (a in 1:min(max_step, i - 1)) for (b in 1:min(max_step, j - 1)) {
      k <- i - a; l <- j - b
      if (!is.finite(E[k, l])) next
      ti <- t[k:i]; slope <- (t[j] - t[l]) / (t[i] - t[k])
      q2w <- .lin_interp(t, q2, t[l] + slope * (ti - t[k])) * sqrt(slope)
      d <- (q1[k:i] - q2w)^2
      c0 <- E[k, l] + sum((d[-1] + d[-length(d)]) / 2 * diff(ti))
      if (c0 < best) { best <- c0; bk <- k; bl <- l }
    }
    E[i, j] <- best; Pi[i, j] <- bk; Pj[i, j] <- bl
  }
  pi_ <- N; pj_ <- N; i <- N; j <- N
  while (i > 1 || j > 1) {
    ni <- Pi[i, j]; nj <- Pj[i, j]
    pi_ <- c(ni, pi_); pj_ <- c(nj, pj_); i <- ni; j <- nj
  }
  g <- .lin_interp(t[pi_], t[pj_], t)
  cummax((g - g[1]) / (g[N] - g[1]))            # normalize + enforce monotone
}

.fisher_rao_phase <- function(gamma, t) {
  gp <- diff(gamma) / diff(t); psi <- sqrt(pmax(0, gp))
  inner <- sum(psi * mean(diff(t)))             # <sqrt(gamma'), 1>
  acos(max(-1, min(1, inner)))
}

#' Elastic (SRVF) alignment of two curves
#'
#' Aligns curve \code{y} to curve \code{x} by an optimal time warping computed in
#' the Square-Root Velocity Function representation, separating phase (the
#' warping) from amplitude (the residual shape difference). A pure-R
#' dynamic-programming solver is used by default; if the \pkg{fdasrvf} package is
#' available it is used for higher accuracy.
#'
#' @param x,y Numeric curves of the same length.
#' @param t Optional grid (defaults to an evenly spaced grid on 0..1).
#' @param lambda Non-negative warping-roughness penalty (passed to \pkg{fdasrvf}
#'   when used; the pure-R solver ignores it beyond boundary/monotonicity).
#' @param method Currently \code{"srvf"}.
#' @param max_step Maximum DP step (bounds the local warping slope) for the
#'   pure-R solver (default: 6).
#' @param smooth Smooth the pure-R warping with a monotone spline (default:
#'   \code{FALSE}).
#' @param use_fdasrvf \code{NA} (default) uses \pkg{fdasrvf} if installed;
#'   \code{TRUE} requires it; \code{FALSE} forces the pure-R solver.
#' @return A list with \code{gamma} (the warping), \code{aligned} (\code{y}
#'   warped to \code{x}), \code{amplitude_distance}, \code{phase_distance}, the
#'   SRVFs, and the solver used.
#' @references
#' Srivastava, A., Wu, W., Kurtek, S., Klassen, E., & Marron, J. S. (2011).
#' Registration of functional data using Fisher-Rao metric. arXiv:1103.3817.
#' @seealso [srvfMean()], [warpApply()]
#' @export
#' @examples
#' t <- seq(0, 1, length.out = 80)
#' x <- exp(-((t - 0.5)^2) / (2 * 0.08^2))
#' y <- exp(-((t - 0.65)^2) / (2 * 0.08^2))
#' al <- elasticAlign(x, y)
#' al$amplitude_distance
elasticAlign <- function(x, y, t = NULL, lambda = 0, method = c("srvf"),
                         max_step = 6L, smooth = FALSE, use_fdasrvf = NA) {
  method <- match.arg(method)
  stopifnot(length(x) == length(y))
  n <- length(x)
  if (is.null(t)) t <- seq(0, 1, length.out = n)

  use_fd <- if (is.na(use_fdasrvf)) requireNamespace("fdasrvf", quietly = TRUE)
            else isTRUE(use_fdasrvf)
  if (isTRUE(use_fdasrvf) && !requireNamespace("fdasrvf", quietly = TRUE)) {
    stop("use_fdasrvf = TRUE but the 'fdasrvf' package is not installed.",
         call. = FALSE)
  }

  if (use_fd) {
    fit <- fdasrvf::pair_align_functions(x, y, t, lambda = lambda)
    gamma <- as.numeric(fit$gam)
    aligned <- as.numeric(fit$f2n)
    solver <- "fdasrvf"
  } else {
    q1 <- .srvf(x, t); q2 <- .srvf(y, t)
    gamma <- .dp_warp(q1, q2, t, max_step)
    if (smooth) {
      sm <- stats::smooth.spline(t, gamma, spar = 0.6)$y
      gamma <- cummax((sm - sm[1]) / (sm[length(sm)] - sm[1]))
    }
    aligned <- .lin_interp(t, y, gamma)
    solver <- "dp"
  }
  q1 <- .srvf(x, t); q2a <- .warp_srvf(.srvf(y, t), gamma, t)
  list(gamma = gamma, aligned = aligned,
       amplitude_distance = sqrt(sum((q1 - q2a)^2) * mean(diff(t))),
       phase_distance = .fisher_rao_phase(gamma, t),
       srvf_x = q1, srvf_y_aligned = q2a, method = method, solver = solver)
}

#' Elastic (Karcher) mean of a set of curves
#'
#' Computes the elastic template mean of repeated curves by iteratively aligning
#' every curve to the current mean and averaging the aligned curves. Because the
#' phase variability is removed before averaging, peaks are preserved rather than
#' blurred out (unlike a naive cross-curve mean).
#'
#' @param curves A list of numeric curves of the same length, or a matrix whose
#'   columns are curves.
#' @param t Optional grid (defaults to 0..1).
#' @param iterations Number of refinement iterations (default: 5).
#' @param max_step Maximum DP step for the pure-R solver (default: 6).
#' @param use_fdasrvf As in [elasticAlign()].
#' @return A list with the elastic \code{mean} curve, the \code{warpings} (a
#'   matrix, one column per curve), and the \code{aligned} curves.
#' @seealso [elasticAlign()]
#' @export
#' @examples
#' t <- seq(0, 1, length.out = 80)
#' curves <- lapply(c(-0.1, 0, 0.1), function(s)
#'   exp(-((t - (0.5 + s))^2) / (2 * 0.05^2)))
#' m <- srvfMean(curves)
#' max(m$mean)
srvfMean <- function(curves, t = NULL, iterations = 5L, max_step = 6L,
                     use_fdasrvf = NA) {
  if (is.matrix(curves)) curves <- lapply(seq_len(ncol(curves)), function(j) curves[, j])
  M <- length(curves); n <- length(curves[[1]])
  if (is.null(t)) t <- seq(0, 1, length.out = n)
  mu <- curves[[1]]
  aligned <- curves
  warps <- matrix(t, n, M)
  for (it in seq_len(iterations)) {
    res <- lapply(curves, function(f)
      elasticAlign(mu, f, t, max_step = max_step, use_fdasrvf = use_fdasrvf))
    aligned <- lapply(res, `[[`, "aligned")
    warps <- vapply(res, `[[`, numeric(n), "gamma")
    mu <- Reduce(`+`, aligned) / M
  }
  list(mean = mu, warpings = warps, aligned = aligned)
}

#' Apply a computed warping to a companion signal
#'
#' Carries a warping computed by [elasticAlign()] onto another channel or
#' modality (e.g. warp an EMG envelope by the diffeomorphism that aligns the
#' EEG).
#'
#' @param signal Numeric signal to warp (same length as the warping grid).
#' @param gamma A warping from [elasticAlign()].
#' @param t Optional grid (defaults to 0..1).
#' @return The warped signal.
#' @seealso [elasticAlign()]
#' @export
#' @examples
#' t <- seq(0, 1, length.out = 50)
#' warpApply(sin(2 * pi * t), gamma = t^1.5, t = t)[1:3]
warpApply <- function(signal, gamma, t = NULL) {
  if (is.null(t)) t <- seq(0, 1, length.out = length(signal))
  .lin_interp(t, signal, gamma)
}

# Elastically warp each PhysioExperiment (2nd onward) to the first using its
# first channel; used by alignSignals(method = "elastic"). Only warps objects
# whose length already matches the reference (same duration after rate
# alignment); the warping is stored in metadata()$elastic_warp.
.warp_pes_to_reference <- function(pes) {
  ref_curve <- SummarizedExperiment::assay(pes[[1]], 1L)[, 1]
  n <- length(ref_curve); t <- seq(0, 1, length.out = n)
  for (i in seq_along(pes)[-1]) {
    A <- as.matrix(SummarizedExperiment::assay(pes[[i]], 1L))
    if (nrow(A) != n) next
    al <- elasticAlign(ref_curve, A[, 1], t)
    Aw <- apply(A, 2, function(col) warpApply(col, al$gamma, t))
    SummarizedExperiment::assay(pes[[i]], 1L) <- Aw
    md <- S4Vectors::metadata(pes[[i]]); md$elastic_warp <- al$gamma
    S4Vectors::metadata(pes[[i]]) <- md
  }
  pes
}
