# Whole-block cross-modal association: how strongly two (or many) modalities are
# related as a WHOLE, beyond a single pair of directions.
#
#   * rvCoefficient()          -- the RV coefficient, a multivariate correlation
#                                 between two blocks (Robert & Escoufier 1976).
#   * distanceCorrelation()    -- distance correlation, sensitive to NON-linear
#                                 dependence, zero iff independent (Szekely 2007).
#   * representationalSimilarity() -- RSA: compare the observation-geometry
#                                 (dissimilarity matrix) of each modality.
# Dependency-free base R, with permutation tests.

.assoc_perm_p <- function(stat_fn, obs, permute, n_perm, seed) {
  if (n_perm <= 0) return(NA_real_)
  if (!is.null(seed)) set.seed(seed)
  n <- permute()
  ge <- sum(vapply(seq_len(n_perm), function(b) stat_fn(sample.int(n)) >= obs, logical(1)))
  (ge + 1) / (n_perm + 1)
}

#' RV coefficient between two modality blocks
#'
#' The multivariate generalisation of a squared correlation: how much the two
#' blocks' configurations of observations agree, in \[0, 1] (0 = unrelated, 1 =
#' identical up to rotation/scaling).
#'
#' @param X,Y Matrices with the same `n` rows (observations).
#' @param center Column-centre each block first (default `TRUE`).
#' @param permutations Row permutations for the null p-value (default 0 = skip).
#' @param seed Optional RNG seed for the permutation test.
#' @return an `rv_coefficient` object: `rv`, `p_value`, `n`.
#' @references Robert P, Escoufier Y (1976) Appl Stat 25:257-265.
#' @seealso [distanceCorrelation()], [representationalSimilarity()], [cca()]
#' @export
#' @examples
#' set.seed(1); z <- matrix(rnorm(40 * 2), 40, 2)
#' X <- z %*% matrix(rnorm(2 * 5), 2, 5) + matrix(rnorm(40 * 5, 0, .3), 40, 5)
#' Y <- z %*% matrix(rnorm(2 * 4), 2, 4) + matrix(rnorm(40 * 4, 0, .3), 40, 4)
#' rvCoefficient(X, Y)$rv
rvCoefficient <- function(X, Y, center = TRUE, permutations = 0L, seed = NULL) {
  X <- as.matrix(X); Y <- as.matrix(Y); n <- nrow(X)
  if (nrow(Y) != n) stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  if (center) { X <- scale(X, center = TRUE, scale = FALSE)
                Y <- scale(Y, center = TRUE, scale = FALSE) }
  rv_of <- function(Yp) {
    sum(crossprod(X, Yp)^2) /
      sqrt(sum(crossprod(X)^2) * sum(crossprod(Yp)^2))
  }
  rv <- rv_of(Y)
  p <- .assoc_perm_p(function(ord) rv_of(Y[ord, , drop = FALSE]), rv,
                     function() n, permutations, seed)
  structure(list(rv = rv, p_value = p, n = n), class = "rv_coefficient")
}

#' @export
print.rv_coefficient <- function(x, ...) {
  cat(sprintf("RV coefficient = %.3f (n = %d)%s\n", x$rv, x$n,
              if (is.na(x$p_value)) "" else sprintf(", p = %.3f", x$p_value)))
  invisible(x)
}

# double-centred Euclidean distance matrix of the rows of M
.assoc_dcenter <- function(M) {
  D <- as.matrix(stats::dist(M))
  rm <- rowMeans(D); gm <- mean(D)
  D - outer(rm, rm, `+`) + gm
}

#' Distance correlation between two modality blocks
#'
#' A dependence measure that captures NON-linear as well as linear cross-modal
#' relationships and is zero if and only if the modalities are independent
#' (Szekely et al. 2007). In \[0, 1].
#'
#' @param X,Y Matrices with the same `n` rows.
#' @param permutations Row permutations for the null p-value (default 0).
#' @param seed Optional RNG seed.
#' @return a `distance_correlation` object: `dcor`, `dcov`, `p_value`, `n`.
#' @references Szekely GJ, Rizzo ML, Bakirov NK (2007) Ann Stat 35:2769-2794.
#' @seealso [rvCoefficient()], [representationalSimilarity()]
#' @export
#' @examples
#' set.seed(1); x <- rnorm(60); y <- x^2 + rnorm(60, 0, 0.3)   # nonlinear
#' distanceCorrelation(x, y)$dcor
distanceCorrelation <- function(X, Y, permutations = 0L, seed = NULL) {
  X <- as.matrix(X); Y <- as.matrix(Y); n <- nrow(X)
  if (nrow(Y) != n) stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  A <- .assoc_dcenter(X); B <- .assoc_dcenter(Y)
  dvarX <- sqrt(mean(A * A)); dvarY <- sqrt(mean(B * B))
  dcov_of <- function(Bp) sqrt(max(mean(A * Bp), 0))
  dcov <- dcov_of(B)
  denom <- sqrt(dvarX * dvarY)
  dcor <- if (denom > 0) dcov / denom else 0
  p <- .assoc_perm_p(function(ord) {
    Bp <- B[ord, ord, drop = FALSE]; dcov_of(Bp) / denom
  }, dcor, function() n, permutations, seed)
  structure(list(dcor = min(dcor, 1), dcov = dcov, p_value = p, n = n),
            class = "distance_correlation")
}

#' @export
print.distance_correlation <- function(x, ...) {
  cat(sprintf("Distance correlation = %.3f (n = %d)%s\n", x$dcor, x$n,
              if (is.na(x$p_value)) "" else sprintf(", p = %.3f", x$p_value)))
  invisible(x)
}

#' Representational similarity across modalities (RSA)
#'
#' Compares the *geometry of observations* between modalities: each block's
#' representational dissimilarity matrix (RDM -- pairwise distances between the
#' shared observations, using that modality's features) is compared across blocks.
#' Two modalities are "representationally similar" if observations that are far
#' apart in one are far apart in the other.
#'
#' @param blocks A list of >= 2 matrices, each `n x p_k` (shared rows).
#' @param method RDM distance: `"euclidean"` (default) or `"correlation"`
#'   (`1 - cor`).
#' @param compare RDM comparison: `"spearman"` (default) or `"pearson"`.
#' @param permutations For a PAIR of blocks, Mantel permutations for a p-value
#'   (default 0).
#' @param seed Optional RNG seed.
#' @return an `rsa` object: `similarity` (block x block RSA matrix) and, for two
#'   blocks with `permutations > 0`, `p_value` (Mantel).
#' @references Kriegeskorte N, et al. (2008) Front Syst Neurosci 2:4.
#' @seealso [rvCoefficient()], [distanceCorrelation()]
#' @export
#' @examples
#' set.seed(1); z <- matrix(rnorm(30 * 2), 30, 2)
#' a <- z %*% matrix(rnorm(2 * 5), 2, 5); b <- z %*% matrix(rnorm(2 * 4), 2, 4)
#' representationalSimilarity(list(a = a, b = b))$similarity
representationalSimilarity <- function(blocks, method = c("euclidean", "correlation"),
                                       compare = c("spearman", "pearson"),
                                       permutations = 0L, seed = NULL) {
  method <- match.arg(method); compare <- match.arg(compare)
  if (!is.list(blocks) || length(blocks) < 2L)
    stop("`blocks` must be a list of >= 2 matrices sharing rows.", call. = FALSE)
  n <- nrow(blocks[[1]])
  if (length(unique(vapply(blocks, nrow, integer(1)))) != 1L)
    stop("all blocks must share rows (observations).", call. = FALSE)
  rdm <- function(M) {
    if (method == "euclidean") as.matrix(stats::dist(M))
    else 1 - stats::cor(t(M))
  }
  low <- function(D) D[lower.tri(D)]
  R <- lapply(blocks, rdm); v <- lapply(R, low)
  K <- length(blocks)
  S <- matrix(1, K, K, dimnames = list(names(blocks), names(blocks)))
  for (i in seq_len(K)) for (j in seq_len(K)) if (i < j) {
    s <- stats::cor(v[[i]], v[[j]], method = compare); S[i, j] <- S[j, i] <- s
  }
  pval <- NA_real_
  if (K == 2L && permutations > 0L) {
    if (!is.null(seed)) set.seed(seed)
    obs <- S[1, 2]; D2 <- R[[2]]
    ge <- sum(vapply(seq_len(permutations), function(b) {
      o <- sample.int(n); Dp <- D2[o, o]
      stats::cor(v[[1]], low(Dp), method = compare) >= obs
    }, logical(1)))
    pval <- (ge + 1) / (permutations + 1)
  }
  structure(list(similarity = S, p_value = pval, method = method,
                 compare = compare, n = n), class = "rsa")
}

#' @export
print.rsa <- function(x, ...) {
  cat(sprintf("Representational similarity (%s RDM, %s) -- %d observations\n",
              x$method, x$compare, x$n))
  print(round(x$similarity, 3))
  if (!is.na(x$p_value)) cat(sprintf("  Mantel p = %.3f\n", x$p_value))
  invisible(x)
}
