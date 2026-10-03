# Two-block cross-modal relationships: which patterns in one modality co-vary
# with which patterns in another (e.g. muscle activations <-> joint kinematics).
#
#   * cca()       -- canonical correlation analysis (regularised), the pairs of
#                    directions with maximal cross-modal correlation.
#   * plsBlocks() -- two-block partial least squares (PLS-SVD), the pairs of
#                    directions with maximal cross-modal COVARIANCE (robust when
#                    features outnumber observations).
# Dependency-free base R.

# symmetric inverse square root with ridge regularisation
.tb_inv_sqrt <- function(C, lambda) {
  e <- eigen(C + diag(lambda, nrow(C)), symmetric = TRUE)
  d <- e$values; d[d < 1e-10] <- 1e-10
  e$vectors %*% (t(e$vectors) / sqrt(d))
}
# `%||%` is defined in fusion-multiblock.R (package-internal)

#' Canonical correlation analysis between two modalities
#'
#' Finds the linear combinations of two modality blocks whose scores are maximally
#' correlated -- the directions in which the modalities co-vary. Ridge
#' regularisation keeps it well-posed when a block has more features than
#' observations.
#'
#' @param X,Y Matrices with the same `n` rows (observations); columns are the
#'   features of each modality.
#' @param ncomp Number of canonical pairs (default: `min(ncol(X), ncol(Y))`).
#' @param lambda_x,lambda_y Ridge penalties; if `0` (default) a small ridge is
#'   added automatically when a block is rank-deficient (features >= n).
#' @return a `cca` object: `cor` (canonical correlations, decreasing), `xcoef`,
#'   `ycoef` (canonical vectors), `xscores`, `yscores` (canonical variates).
#' @references Hotelling H (1936); Vinod HD (1976) ridge CCA.
#' @seealso [plsBlocks()], [rvCoefficient()]
#' @export
#' @examples
#' set.seed(1); z <- matrix(rnorm(50 * 2), 50, 2)
#' X <- z %*% matrix(rnorm(2 * 5), 2, 5) + matrix(rnorm(50 * 5, 0, .3), 50, 5)
#' Y <- z %*% matrix(rnorm(2 * 4), 2, 4) + matrix(rnorm(50 * 4, 0, .3), 50, 4)
#' cca(X, Y)$cor
cca <- function(X, Y, ncomp = NULL, lambda_x = 0, lambda_y = 0) {
  X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
  Y <- scale(as.matrix(Y), center = TRUE, scale = FALSE)
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  if (nrow(Y) != n) stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  Cxx <- crossprod(X) / (n - 1); Cyy <- crossprod(Y) / (n - 1)
  Cxy <- crossprod(X, Y) / (n - 1)
  if (lambda_x == 0 && p >= n) lambda_x <- 1e-3 * mean(diag(Cxx))   # auto-ridge
  if (lambda_y == 0 && q >= n) lambda_y <- 1e-3 * mean(diag(Cyy))
  ix <- .tb_inv_sqrt(Cxx, lambda_x); iy <- .tb_inv_sqrt(Cyy, lambda_y)
  s <- svd(ix %*% Cxy %*% iy)
  nc <- min(ncomp %||% length(s$d), length(s$d), p, q)
  A <- ix %*% s$u[, seq_len(nc), drop = FALSE]
  B <- iy %*% s$v[, seq_len(nc), drop = FALSE]
  structure(list(cor = pmin(s$d[seq_len(nc)], 1), xcoef = A, ycoef = B,
                 xscores = X %*% A, yscores = Y %*% B, n = n, ncomp = nc),
            class = "cca")
}

#' @export
print.cca <- function(x, ...) {
  cat(sprintf("Canonical correlation -- %d pairs, %d observations\n", x$ncomp, x$n))
  cat("  canonical correlations:", paste(sprintf("%.3f", x$cor), collapse = " "), "\n")
  invisible(x)
}

#' Two-block partial least squares (PLS-SVD)
#'
#' Finds the linear combinations of two modality blocks with maximal cross-modal
#' COVARIANCE (unlike CCA's correlation), by the singular value decomposition of
#' the cross-covariance and successive deflation. Robust when features outnumber
#' observations, so it suits high-dimensional modalities.
#'
#' @param X,Y Matrices with the same `n` rows (observations).
#' @param ncomp Number of PLS components (default 2).
#' @return a `pls_blocks` object: `xweights`, `yweights` (the paired directions),
#'   `xscores`, `yscores`, `covariance` (per-component score covariance) and
#'   `correlation` (per-component score correlation).
#' @references Wold H (1975); Bookstein (1994) PLS-SVD; Krishnan et al. (2011).
#' @seealso [cca()], [coupledNMF()]
#' @export
#' @examples
#' set.seed(1); z <- matrix(rnorm(40), 40, 1)
#' X <- z %*% matrix(rnorm(6), 1, 6) + matrix(rnorm(40 * 6, 0, .3), 40, 6)
#' Y <- z %*% matrix(rnorm(5), 1, 5) + matrix(rnorm(40 * 5, 0, .3), 40, 5)
#' plsBlocks(X, Y, ncomp = 1)$correlation
plsBlocks <- function(X, Y, ncomp = 2L) {
  X0 <- scale(as.matrix(X), center = TRUE, scale = FALSE)
  Y0 <- scale(as.matrix(Y), center = TRUE, scale = FALSE)
  n <- nrow(X0); if (nrow(Y0) != n) stop("`X` and `Y` must have the same rows.", call. = FALSE)
  nc <- min(ncomp, ncol(X0), ncol(Y0))
  Xd <- X0; Yd <- Y0
  W <- matrix(0, ncol(X0), nc); C <- matrix(0, ncol(Y0), nc)
  Tt <- matrix(0, n, nc); U <- matrix(0, n, nc); cov_c <- numeric(nc)
  for (k in seq_len(nc)) {
    s <- svd(crossprod(Xd, Yd), nu = 1, nv = 1)       # first singular vectors of X'Y
    w <- s$u[, 1]; c <- s$v[, 1]
    t <- Xd %*% w; u <- Yd %*% c
    W[, k] <- w; C[, k] <- c; Tt[, k] <- t; U[, k] <- u; cov_c[k] <- s$d[1] / (n - 1)
    Xd <- Xd - t %*% (crossprod(t, Xd) / sum(t^2))    # deflate each block by its score
    Yd <- Yd - u %*% (crossprod(u, Yd) / sum(u^2))
  }
  corr <- vapply(seq_len(nc), function(k) stats::cor(Tt[, k], U[, k]), numeric(1))
  structure(list(xweights = W, yweights = C, xscores = Tt, yscores = U,
                 covariance = cov_c, correlation = corr, n = n, ncomp = nc),
            class = "pls_blocks")
}

#' @export
print.pls_blocks <- function(x, ...) {
  cat(sprintf("Two-block PLS -- %d components, %d observations\n", x$ncomp, x$n))
  cat("  score correlation:", paste(sprintf("%.3f", x$correlation), collapse = " "), "\n")
  invisible(x)
}
