# Multi-block latent fusion of several modality blocks.
#
# The coupling functions relate TWO signals; these integrate MANY modality blocks
# (kinematics, kinetics, EMG, IMU, ...) sharing a common observation dimension
# (e.g. subjects, or gait-cycle points) into a JOINT structure shared across all
# modalities plus INDIVIDUAL structure specific to each -- the integrated
# multi-modal representation of gait. Two established methods, dependency-free:
#   * JIVE (Joint and Individual Variation Explained; Lock et al. 2013) --
#     separates the variation common to all blocks from block-specific variation.
#   * Multiple Factor Analysis (MFA; Escofier & Pages 1994) -- a multi-block PCA
#     that weights each block so none dominates, giving global factors with
#     per-block contributions.

# rank-r reconstruction from an svd object
.mb_lowrank <- function(s, r) {
  r <- min(r, length(s$d)); if (r < 1) return(matrix(0, nrow(s$u), nrow(s$v)))
  s$u[, seq_len(r), drop = FALSE] %*% diag(s$d[seq_len(r)], r) %*% t(s$v[, seq_len(r), drop = FALSE])
}

# coerce a list of blocks (each n x p_k, shared rows) -> list of p_k x n
.mb_blocks <- function(blocks) {
  if (!is.list(blocks) || length(blocks) < 2L)
    stop("`blocks` must be a list of >= 2 matrices sharing rows (observations).", call. = FALSE)
  Xs <- lapply(blocks, function(b) t(as.matrix(b)))
  n <- vapply(Xs, ncol, integer(1))
  if (length(unique(n)) != 1L)
    stop("all blocks must have the same number of rows (shared observations).", call. = FALSE)
  if (is.null(names(Xs))) names(Xs) <- paste0("block", seq_along(Xs))
  Xs
}

#' JIVE: Joint and Individual Variation Explained
#'
#' Decomposes several modality blocks that share an observation dimension into a
#' low-rank JOINT structure common to all blocks, low-rank INDIVIDUAL structure
#' specific to each, and residual noise (Lock et al. 2013). The joint scores are
#' the integrated multi-modal representation; the per-block variance split tells
#' you how much of each modality is shared versus idiosyncratic.
#'
#' @param blocks A list of >= 2 matrices, each `n x p_k`: the same `n`
#'   observations (e.g. subjects, or gait-cycle points), different features per
#'   modality.
#' @param rank_joint Rank of the joint structure.
#' @param rank_individual Rank of each block's individual structure (a scalar,
#'   recycled, or a per-block vector).
#' @param center,scale Row-center each feature and scale each block to unit
#'   Frobenius norm (so blocks contribute comparably); both default `TRUE`.
#' @param max_iter,tol Iteration budget and convergence tolerance.
#' @return a `jive` object: `joint` and `individual` (lists of `n x p_k`
#'   matrices), `joint_scores` (`n x rank_joint`, the shared component scores),
#'   `variance` (data frame of joint / individual / residual fraction per block),
#'   `ranks`, `iterations`.
#' @references Lock EF, Hoadley KA, Marron JS, Nobel AB (2013) Ann Appl Stat
#'   7:523-542.
#' @seealso [multipleFactorAnalysis()], [cca()]
#' @export
#' @examples
#' set.seed(1)
#' shared <- matrix(rnorm(40 * 2), 40, 2)                 # common subject scores
#' k <- shared %*% matrix(rnorm(2 * 6), 2, 6) + matrix(rnorm(40 * 6, 0, .3), 40, 6)
#' e <- shared %*% matrix(rnorm(2 * 5), 2, 5) + matrix(rnorm(40 * 5, 0, .3), 40, 5)
#' fit <- jive(list(kin = k, emg = e), rank_joint = 2, rank_individual = 1)
#' fit$variance
jive <- function(blocks, rank_joint = 1L, rank_individual = 1L,
                 center = TRUE, scale = TRUE, max_iter = 1000L, tol = 1e-8) {
  Xs <- .mb_blocks(blocks); K <- length(Xs); n <- ncol(Xs[[1]])
  if (center) Xs <- lapply(Xs, function(X) X - rowMeans(X))
  fro <- vapply(Xs, function(X) sqrt(sum(X^2)), numeric(1)); fro[fro == 0] <- 1
  Xw <- if (scale) Map(function(X, f) X / f, Xs, fro) else Xs
  pk <- vapply(Xw, nrow, integer(1)); ends <- cumsum(pk); starts <- ends - pk + 1L
  ri <- if (length(rank_individual) == 1L) rep(rank_individual, K) else rank_individual
  A <- lapply(Xw, function(X) matrix(0, nrow(X), ncol(X)))
  Jfull <- matrix(0, sum(pk), n); prev <- Inf
  it <- 0L
  repeat {
    it <- it + 1L
    R <- do.call(rbind, Map(`-`, Xw, A))            # remove individual, estimate joint
    sj <- svd(R); rj <- min(rank_joint, length(sj$d))
    Vj <- sj$v[, seq_len(rj), drop = FALSE]         # n x rj shared sample space
    Jfull <- .mb_lowrank(sj, rj)
    Jk <- lapply(seq_len(K), function(k) Jfull[starts[k]:ends[k], , drop = FALSE])
    Pperp <- diag(n) - tcrossprod(Vj)               # orthocomplement of joint
    A <- lapply(seq_len(K), function(k) {
      M <- (Xw[[k]] - Jk[[k]]) %*% Pperp
      .mb_lowrank(svd(M), ri[k])
    })
    delta <- sum((Jfull - prev)^2) / max(sum(Jfull^2), .Machine$double.eps)
    if (is.finite(delta) && delta < tol) break
    if (it >= max_iter) break
    prev <- Jfull
  }
  var_tab <- data.frame(
    block = names(Xw),
    joint = vapply(seq_len(K), function(k) sum(Jk[[k]]^2) / sum(Xw[[k]]^2), numeric(1)),
    individual = vapply(seq_len(K), function(k) sum(A[[k]]^2) / sum(Xw[[k]]^2), numeric(1)),
    row.names = NULL)
  var_tab$residual <- 1 - var_tab$joint - var_tab$individual
  structure(list(
    joint = stats::setNames(lapply(Jk, t), names(Xw)),
    individual = stats::setNames(lapply(A, t), names(Xw)),
    joint_scores = Vj, variance = var_tab,
    ranks = list(joint = rank_joint, individual = ri), iterations = it,
    n = n, K = K), class = "jive")
}

#' @export
print.jive <- function(x, ...) {
  cat(sprintf("JIVE -- %d blocks, %d shared observations, joint rank %d (%d iters)\n",
              x$K, x$n, x$ranks$joint, x$iterations))
  v <- x$variance
  for (i in seq_len(nrow(v)))
    cat(sprintf("  %-10s joint %4.1f%%  individual %4.1f%%  residual %4.1f%%\n",
                v$block[i], 100 * v$joint[i], 100 * v$individual[i], 100 * v$residual[i]))
  invisible(x)
}

#' Multiple Factor Analysis (multi-block PCA)
#'
#' A global PCA over several blocks in which each block is first scaled by its
#' leading singular value, so a large or high-variance modality cannot dominate
#' the shared factors (Escofier & Pages 1994). Gives global factor scores on the
#' shared observations, with each block's contribution to each factor.
#'
#' @param blocks A list of >= 2 matrices, each `n x p_k` (shared rows).
#' @param ncomp Number of global factors (default 2).
#' @param center,scale_columns Center (and, if `TRUE`, unit-variance scale) each
#'   column before the block weighting.
#' @return an `mfa` object: `scores` (`n x ncomp` global factor scores),
#'   `eigenvalues`, `explained` (proportion of variance per factor),
#'   `contributions` (block x factor, each block's share of a factor's inertia),
#'   `block_weights`.
#' @references Escofier B, Pages J (1994) Comput Stat Data Anal 18:121-140.
#' @seealso [jive()], [rvCoefficient()]
#' @export
#' @examples
#' set.seed(1)
#' s <- matrix(rnorm(30 * 2), 30, 2)
#' b1 <- s %*% matrix(rnorm(2 * 5), 2, 5) + matrix(rnorm(30 * 5, 0, .3), 30, 5)
#' b2 <- s %*% matrix(rnorm(2 * 4), 2, 4) + matrix(rnorm(30 * 4, 0, .3), 30, 4)
#' multipleFactorAnalysis(list(kin = b1, grf = b2))$explained
multipleFactorAnalysis <- function(blocks, ncomp = 2L, center = TRUE,
                                   scale_columns = FALSE) {
  if (!is.list(blocks) || length(blocks) < 2L)
    stop("`blocks` must be a list of >= 2 matrices sharing rows.", call. = FALSE)
  Bs <- lapply(blocks, as.matrix); n <- nrow(Bs[[1]])
  if (length(unique(vapply(Bs, nrow, integer(1)))) != 1L)
    stop("all blocks must share rows (observations).", call. = FALSE)
  Bs <- lapply(Bs, function(B) scale(B, center = center, scale = scale_columns))
  w <- vapply(Bs, function(B) 1 / (svd(B)$d[1]^2), numeric(1))      # 1 / first eigenvalue
  Bw <- Map(function(B, wi) B * sqrt(wi), Bs, w)
  X <- do.call(cbind, Bw)
  sv <- svd(X); nc <- min(ncomp, length(sv$d))
  scores <- sv$u[, seq_len(nc), drop = FALSE] %*% diag(sv$d[seq_len(nc)], nc)
  ev <- sv$d^2; explained <- ev / sum(ev)
  pk <- vapply(Bw, ncol, integer(1)); ends <- cumsum(pk); starts <- ends - pk + 1L
  V <- sv$v[, seq_len(nc), drop = FALSE]
  contrib <- t(vapply(seq_along(Bw), function(k)
    colSums(V[starts[k]:ends[k], , drop = FALSE]^2), numeric(nc)))   # block x factor
  dimnames(contrib) <- list(names(blocks) %||% paste0("block", seq_along(Bw)),
                            paste0("F", seq_len(nc)))
  structure(list(scores = scores, eigenvalues = ev[seq_len(nc)],
                 explained = explained[seq_len(nc)], contributions = contrib,
                 block_weights = w, n = n, ncomp = nc), class = "mfa")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.mfa <- function(x, ...) {
  cat(sprintf("Multiple Factor Analysis -- %d blocks, %d observations, %d factors\n",
              nrow(x$contributions), x$n, x$ncomp))
  cat("  variance explained:", paste(sprintf("%.1f%%", 100 * x$explained), collapse = " "), "\n")
  invisible(x)
}
