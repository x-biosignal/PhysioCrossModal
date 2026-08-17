# Coupled (joint) non-negative factorization of two modality blocks.
#
# When two non-negative modalities are driven by a COMMON latent activation --
# e.g. EMG envelopes and a kinematic-speed envelope sharing the same temporal
# drive -- a coupled factorization forces one shared factor matrix while letting
# each modality keep its own loadings, so the shared drive is estimated from both
# at once. Multiplicative non-negative updates; dependency-free base R.

#' Coupled non-negative matrix factorization across two modalities
#'
#' Factorises two non-negative blocks that share the observation dimension into a
#' SHARED factor matrix `W` and modality-specific loadings `Hx`, `Hy`, minimising
#' `||X - W Hx||^2 + ||Y - W Hy||^2`. The shared `W` is the common latent drive
#' recovered jointly from both modalities.
#'
#' @param X,Y Non-negative matrices with the same `n` rows (shared observations,
#'   e.g. time or gait-cycle points); columns are each modality's channels.
#' @param rank Number of shared factors.
#' @param max_iter,tol Iteration budget and convergence tolerance.
#' @param n_restart Random restarts (best kept; default 5).
#' @param seed Optional RNG seed.
#' @return a `coupled_nmf` object: `shared` (`n x rank` shared factor `W`), `Hx`,
#'   `Hy` (`rank x p` and `rank x q` loadings), `vaf_x`, `vaf_y`, `vaf` (combined),
#'   `iterations`.
#' @references Lee & Seung (2001) NMF; Cichocki et al. (2009) coupled/joint NMF.
#' @seealso [plsBlocks()], [jive()]
#' @export
#' @examples
#' set.seed(1); W <- matrix(runif(50 * 2), 50, 2)      # shared drive
#' X <- W %*% matrix(runif(2 * 6), 2, 6); Y <- W %*% matrix(runif(2 * 4), 2, 4)
#' fit <- coupledNMF(X, Y, rank = 2, n_restart = 2)
#' fit$vaf
coupledNMF <- function(X, Y, rank, max_iter = 500L, tol = 1e-7,
                       n_restart = 5L, seed = NULL) {
  X <- as.matrix(X); Y <- as.matrix(Y); n <- nrow(X)
  if (nrow(Y) != n) stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  if (any(X < 0) || any(Y < 0)) stop("`X` and `Y` must be non-negative.", call. = FALSE)
  p <- ncol(X); q <- ncol(Y); eps <- 1e-10
  normX <- sum(X^2); normY <- sum(Y^2)
  if (!is.null(seed)) set.seed(seed)
  best <- NULL
  for (rs in seq_len(n_restart)) {
    W <- matrix(stats::runif(n * rank), n, rank)
    Hx <- matrix(stats::runif(rank * p), rank, p)
    Hy <- matrix(stats::runif(rank * q), rank, q)
    prev <- Inf
    for (it in seq_len(max_iter)) {
      W  <- W * (X %*% t(Hx) + Y %*% t(Hy)) /
              (W %*% (tcrossprod(Hx) + tcrossprod(Hy)) + eps)
      WtW <- crossprod(W)
      Hx <- Hx * (crossprod(W, X)) / (WtW %*% Hx + eps)
      Hy <- Hy * (crossprod(W, Y)) / (WtW %*% Hy + eps)
      if (it %% 25L == 0L || it == max_iter) {
        resid <- sum((X - W %*% Hx)^2) + sum((Y - W %*% Hy)^2)
        if (is.finite(prev) && abs(prev - resid) < tol * (prev + eps)) break
        prev <- resid
      }
    }
    rx <- sum((X - W %*% Hx)^2); ry <- sum((Y - W %*% Hy)^2)
    if (is.null(best) || (rx + ry) < best$resid)
      best <- list(W = W, Hx = Hx, Hy = Hy, rx = rx, ry = ry, resid = rx + ry, iter = it)
  }
  structure(list(shared = best$W, Hx = best$Hx, Hy = best$Hy,
                 vaf_x = 1 - best$rx / normX, vaf_y = 1 - best$ry / normY,
                 vaf = 1 - best$resid / (normX + normY),
                 rank = rank, iterations = best$iter, n = n), class = "coupled_nmf")
}

#' @export
print.coupled_nmf <- function(x, ...) {
  cat(sprintf("Coupled NMF -- rank %d, %d shared observations (%d iters)\n",
              x$rank, x$n, x$iterations))
  cat(sprintf("  VAF: block X %.3f, block Y %.3f, combined %.3f\n",
              x$vaf_x, x$vaf_y, x$vaf))
  invisible(x)
}
