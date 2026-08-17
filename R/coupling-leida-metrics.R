# LEiDA state-transition metrics ------------------------------------------
# Summarise the dynamics of the LEiDA state sequence: the Markov transition
# matrix between connectivity states, how often the brain switches, how long it
# dwells, and how (un)predictable the switching is.

#' LEiDA state-transition metrics
#'
#' From a [leidaStates()] result, computes the state-to-state transition
#' probability matrix, the switching rate, fractional occupancy, mean dwell
#' time, and the occupancy-weighted transition entropy (predictability of the
#' next state).
#'
#' @param x A `"leida"` object from [leidaStates()].
#' @param sr Optional sampling rate (Hz); if given, dwell times are returned in
#'   seconds and the switching rate additionally in Hz.
#' @return An object of class `"leida_transitions"`: a list with `transition`
#'   (row-stochastic `state x state` matrix), `counts` (raw transition counts),
#'   `switching_rate` (fraction of steps that change state), `switching_rate_hz`
#'   (switches per second, or `NA`), `occupancy`, `dwell`, and `entropy`
#'   (mean bits of uncertainty about the next state).
#' @examples
#' set.seed(1); sr <- 100; n <- 2000; t <- seq_len(n) / sr
#' base <- sin(2 * pi * 10 * t)
#' X <- cbind(base + rnorm(n, sd = .3), base + rnorm(n, sd = .3),
#'            sin(2 * pi * 10 * t + 1) + rnorm(n, sd = .3),
#'            sin(2 * pi * 10 * t + 1) + rnorm(n, sd = .3))
#' res <- leidaStates(X, sr = sr, freq_band = c(8, 12), n_states = 2, seed = 1)
#' leidaTransitions(res, sr = sr)
#' @export
leidaTransitions <- function(x, sr = NULL) {
  stopifnot(inherits(x, "leida"))
  s <- x$states; K <- x$n_states
  from <- s[-length(s)]; to <- s[-1]

  counts <- matrix(0L, K, K)
  for (k in seq_along(from)) counts[from[k], to[k]] <- counts[from[k], to[k]] + 1L
  rs <- rowSums(counts); denom <- ifelse(rs == 0, 1, rs)
  transition <- counts / denom

  n_switch <- sum(from != to)
  switching_rate <- if (length(from)) n_switch / length(from) else 0
  switching_rate_hz <- if (!is.null(sr)) n_switch / (length(s) / sr) else NA_real_
  dwell <- if (!is.null(sr)) x$dwell / sr else x$dwell

  ent <- function(p) { p <- p[p > 0]; if (!length(p)) 0 else -sum(p * log2(p)) }
  entropy <- sum(apply(transition, 1L, ent) * x$occupancy)

  structure(
    list(transition = transition, counts = counts,
         switching_rate = switching_rate, switching_rate_hz = switching_rate_hz,
         occupancy = x$occupancy, dwell = dwell, entropy = entropy,
         n_states = K),
    class = "leida_transitions")
}

#' @export
print.leida_transitions <- function(x, ...) {
  cat("<LEiDA state-transition metrics>\n")
  cat(sprintf("  states:         %d\n", x$n_states))
  cat(sprintf("  switching rate: %.3f per step%s\n", x$switching_rate,
              if (!is.na(x$switching_rate_hz))
                sprintf(" (%.2f Hz)", x$switching_rate_hz) else ""))
  cat(sprintf("  transition entropy: %.2f bits\n", x$entropy))
  invisible(x)
}
