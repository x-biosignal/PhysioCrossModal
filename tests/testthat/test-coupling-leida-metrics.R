# build a leida object directly from a known state sequence
fake_leida <- function(states, K = max(states)) {
  r <- rle(states)
  dwell <- vapply(seq_len(K), function(s) {
    runs <- r$lengths[r$values == s]; if (length(runs)) mean(runs) else 0
  }, numeric(1))
  structure(list(states = states, occupancy = tabulate(states, K) / length(states),
                 dwell = dwell, n_states = K,
                 centroids = matrix(0, K, 2)), class = "leida")
}

test_that("transition matrix and switching rate match a known sequence", {
  # deterministic 1 -> 2 -> 1 -> 2 ... : always switches
  lo <- fake_leida(rep(c(1, 2), 50))
  tr <- leidaTransitions(lo)
  expect_s3_class(tr, "leida_transitions")
  expect_equal(tr$transition[1, 2], 1)      # 1 always -> 2
  expect_equal(tr$transition[2, 1], 1)      # 2 always -> 1
  expect_equal(tr$switching_rate, 1)        # every step switches
  expect_equal(tr$entropy, 0)               # deterministic -> 0 bits
})

test_that("dwell and switching reflect long runs", {
  # long runs: 1x20 then 2x20, repeated -> rare switching
  s <- rep(rep(c(1, 2), each = 20), 5)
  tr <- leidaTransitions(fake_leida(s), sr = 10)
  expect_lt(tr$switching_rate, 0.1)
  # dwell ~ 20 samples = 2 s at 10 Hz
  expect_equal(mean(tr$dwell), 2, tolerance = 0.3)
  expect_false(is.na(tr$switching_rate_hz))
})

test_that("transition rows are stochastic and entropy is positive when random", {
  set.seed(1)
  s <- sample(1:3, 900, replace = TRUE)     # random -> near-max entropy
  tr <- leidaTransitions(fake_leida(s, K = 3))
  expect_equal(rowSums(tr$transition), rep(1, 3), tolerance = 1e-8)
  expect_gt(tr$entropy, 1)                   # close to log2(3) ~ 1.58
})

test_that("works end-to-end on a real leidaStates result", {
  set.seed(1); sr <- 100; n <- 2000; t <- seq_len(n) / sr
  base <- sin(2 * pi * 10 * t)
  X <- cbind(base + rnorm(n, sd = .3), base + rnorm(n, sd = .3),
             sin(2 * pi * 10 * t + 1) + rnorm(n, sd = .3),
             sin(2 * pi * 10 * t + 1) + rnorm(n, sd = .3))
  res <- leidaStates(X, sr = sr, freq_band = c(8, 12), n_states = 2, seed = 1)
  tr <- leidaTransitions(res, sr = sr)
  expect_equal(dim(tr$transition), c(2, 2))
  expect_true(tr$switching_rate >= 0 && tr$switching_rate <= 1)
})
