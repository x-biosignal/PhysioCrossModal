# A composite multi-modal gait index, feature-level fusion, and a cross-block
# latent-factor model -- turning several aligned modalities into one number, one
# design matrix, or one shared latent structure.
#
#   * multimodalGaitIndex()  -- a Gait-Deviation-Index-style composite deviation
#                               from a normative reference, but spanning modalities
#                               (kinematics + kinetics + EMG + spatiotemporal),
#                               with a per-modality breakdown.
#   * fuseBlocks()           -- feature-level fusion into a single block-balanced
#                               design matrix for downstream (multi-view) ML.
#   * crossBlockFactor()     -- a common-factor model over the pooled blocks,
#                               flagging which latent factors are shared across
#                               modalities versus modality-specific.
# Dependency-free base R (+ stats::factanal for the factor model).

# GDI-style standardized-PC distance from a reference cloud
.idx_pc_distance <- function(Xc, Rc, ve) {
  sv <- svd(Rc); ev <- sv$d^2
  nc <- max(1L, which(cumsum(ev) / sum(ev) >= ve)[1])
  V <- sv$v[, seq_len(nc), drop = FALSE]
  sdpc <- sv$d[seq_len(nc)] / sqrt(nrow(Rc) - 1); sdpc[sdpc < 1e-9] <- 1e-9
  list(x = sqrt(rowSums(sweep(Xc %*% V, 2L, sdpc, "/")^2)),
       ref = sqrt(rowSums(sweep(Rc %*% V, 2L, sdpc, "/")^2)))
}

#' Composite multi-modal gait deviation index
#'
#' A Gait-Deviation-Index-style score (Schwartz & Rozumalski 2008) generalised
#' across modalities: how far each subject's multi-modal gait feature vector lies
#' from a normative reference cloud, scaled so a typical control scores ~100 and
#' every 10 points is one standard deviation of deviation. Optionally broken down
#' per modality to show which modality drives the deviation.
#'
#' @param features An `n x p` matrix of subjects' multi-modal gait features
#'   (kinematics, kinetics, EMG summaries, spatiotemporal, ... concatenated).
#' @param reference An `m x p` matrix of the same features for a normative
#'   (control) sample.
#' @param blocks Optional named list mapping each modality to its column indices
#'   in `features`, for a per-modality breakdown.
#' @param ve Reference variance retained by the PCA whitening (default 0.95).
#' @param mean_index,sd_index Scaling (default 100 and 10, as for GDI).
#' @return a `multimodal_index` object: `index` (per subject; ~100 = typical),
#'   `distance`, `reference_mean_index`, and (if `blocks` given) `by_modality`
#'   (subject x modality index matrix).
#' @references Schwartz MH, Rozumalski A (2008) Gait Posture 28:351-357.
#' @seealso [jive()], [fuseBlocks()]
#' @export
#' @examples
#' set.seed(1); ctrl <- matrix(rnorm(60 * 8), 60, 8)
#' subj <- rbind(matrix(rnorm(5 * 8), 5, 8),               # typical
#'               matrix(rnorm(5 * 8, 3), 5, 8))            # deviant
#' multimodalGaitIndex(subj, ctrl)$index
multimodalGaitIndex <- function(features, reference, blocks = NULL, ve = 0.95,
                                mean_index = 100, sd_index = 10) {
  X <- as.matrix(features); Rf <- as.matrix(reference)
  if (ncol(X) != ncol(Rf)) stop("`features` and `reference` need the same columns.", call. = FALSE)
  cm <- colMeans(Rf)
  gdi <- function(cols) {
    Rc <- sweep(Rf[, cols, drop = FALSE], 2L, cm[cols])
    Xc <- sweep(X[, cols, drop = FALSE], 2L, cm[cols])
    d <- .idx_pc_distance(Xc, Rc, ve)
    lr <- log(d$ref + 1e-9); m <- mean(lr); s <- stats::sd(lr); if (s < 1e-9) s <- 1
    list(index = mean_index - sd_index * (log(d$x + 1e-9) - m) / s,
         ref_index = mean_index - sd_index * (lr - m) / s, distance = d$x)
  }
  all <- gdi(seq_len(ncol(X)))
  by_mod <- NULL
  if (!is.null(blocks)) {
    by_mod <- vapply(blocks, function(cols) gdi(cols)$index, numeric(nrow(X)))
    if (is.null(dim(by_mod))) by_mod <- matrix(by_mod, nrow = nrow(X))
    colnames(by_mod) <- names(blocks)
  }
  structure(list(index = all$index, distance = all$distance,
                 reference_mean_index = mean(all$ref_index), by_modality = by_mod,
                 n = nrow(X)), class = "multimodal_index")
}

#' @export
print.multimodal_index <- function(x, ...) {
  cat(sprintf("Multi-modal gait deviation index -- %d subjects\n", x$n))
  cat(sprintf("  index range [%.1f, %.1f]  (reference ~ %.0f)\n",
              min(x$index), max(x$index), x$reference_mean_index))
  invisible(x)
}

#' Feature-level fusion of modality blocks
#'
#' Concatenates several modality blocks into a single design matrix for
#' downstream (multi-view) machine learning, scaling each block so a large or
#' high-variance modality does not dominate the fused features.
#'
#' @param blocks A list of matrices, each `n x p_k` (shared rows).
#' @param method Block scaling before concatenation: `"zscore"` (per-column
#'   standardise, default), `"mfa"` (divide each block by its leading singular
#'   value), or `"frobenius"` (unit Frobenius norm per block).
#' @return a `fused_blocks` object: `matrix` (`n x sum(p_k)`) and `block_index`
#'   (column ranges per block).
#' @seealso [multipleFactorAnalysis()], [multimodalGaitIndex()]
#' @export
#' @examples
#' fuseBlocks(list(a = matrix(rnorm(20), 10, 2), b = matrix(rnorm(30) * 50, 10, 3)))$matrix
fuseBlocks <- function(blocks, method = c("zscore", "mfa", "frobenius")) {
  method <- match.arg(method)
  if (!is.list(blocks) || length(blocks) < 2L)
    stop("`blocks` must be a list of >= 2 matrices sharing rows.", call. = FALSE)
  Bs <- lapply(blocks, as.matrix)
  if (length(unique(vapply(Bs, nrow, integer(1)))) != 1L)
    stop("all blocks must share rows.", call. = FALSE)
  Bs <- switch(method,
    zscore = lapply(Bs, function(B) scale(B)),
    mfa = lapply(Bs, function(B) { Bc <- scale(B, scale = FALSE); Bc / svd(Bc)$d[1] }),
    frobenius = lapply(Bs, function(B) { Bc <- scale(B, scale = FALSE); Bc / sqrt(sum(Bc^2)) }))
  pk <- vapply(Bs, ncol, integer(1)); ends <- cumsum(pk); starts <- ends - pk + 1L
  M <- do.call(cbind, Bs)
  bi <- data.frame(block = names(blocks) %||% paste0("block", seq_along(Bs)),
                   start = starts, end = ends, row.names = NULL)
  structure(list(matrix = M, block_index = bi, method = method), class = "fused_blocks")
}

#' @export
print.fused_blocks <- function(x, ...) {
  cat(sprintf("Fused blocks (%s) -- %d x %d, %d modalities\n",
              x$method, nrow(x$matrix), ncol(x$matrix), nrow(x$block_index)))
  invisible(x)
}

#' Cross-block common-factor model
#'
#' Fits a maximum-likelihood common-factor model to the pooled, standardised
#' modality blocks and flags which latent factors are SHARED across modalities
#' (load on more than one block) versus modality-specific -- a latent-variable
#' view of multi-modal integration.
#'
#' @param blocks A list of >= 2 matrices, each `n x p_k` (shared rows). Requires
#'   more observations than total features.
#' @param nfactors Number of common factors.
#' @param load_threshold A factor is considered to load on a block if any of the
#'   block's variables has `|loading| >=` this (default 0.4).
#' @return a `crossblock_factor` object: `loadings` (variables x factors, with a
#'   `block` attribute), `shared` (logical per factor: loads on >= 2 blocks),
#'   `n_shared`, and the fitted `factanal`.
#' @references Spearman/Thurstone common-factor model; McDonald (1985).
#' @seealso [jive()], [multipleFactorAnalysis()]
#' @export
#' @examples
#' set.seed(1); f <- matrix(rnorm(80 * 1), 80, 1)
#' a <- f %*% matrix(rnorm(4), 1, 4) + matrix(rnorm(80 * 4, 0, .4), 80, 4)
#' b <- f %*% matrix(rnorm(3), 1, 3) + matrix(rnorm(80 * 3, 0, .4), 80, 3)
#' crossBlockFactor(list(a = a, b = b), nfactors = 1)$shared
crossBlockFactor <- function(blocks, nfactors, load_threshold = 0.4) {
  if (!is.list(blocks) || length(blocks) < 2L)
    stop("`blocks` must be a list of >= 2 matrices sharing rows.", call. = FALSE)
  Bs <- lapply(blocks, as.matrix); n <- nrow(Bs[[1]])
  if (length(unique(vapply(Bs, nrow, integer(1)))) != 1L)
    stop("all blocks must share rows.", call. = FALSE)
  block_of <- rep(names(blocks) %||% paste0("block", seq_along(Bs)),
                  vapply(Bs, ncol, integer(1)))
  X <- do.call(cbind, lapply(Bs, scale))
  colnames(X) <- make.unique(paste0(block_of, ".", unlist(lapply(Bs, function(B)
    seq_len(ncol(B))))))
  fa <- stats::factanal(X, factors = nfactors, scores = "regression")
  L <- unclass(fa$loadings); attr(L, "block") <- block_of
  # which factors load on >= 2 distinct blocks
  shared <- vapply(seq_len(nfactors), function(j) {
    blk_hit <- unique(block_of[abs(L[, j]) >= load_threshold])
    length(blk_hit) >= 2L
  }, logical(1))
  structure(list(loadings = L, block = block_of, shared = shared,
                 n_shared = sum(shared), nfactors = nfactors, factanal = fa,
                 scores = fa$scores), class = "crossblock_factor")
}

#' @export
print.crossblock_factor <- function(x, ...) {
  cat(sprintf("Cross-block factor model -- %d factors, %d shared across modalities\n",
              x$nfactors, x$n_shared))
  invisible(x)
}
