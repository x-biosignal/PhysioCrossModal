test_that("composite keys recover correspondence and joint subspace", {
  set.seed(230923)
  s <- rnorm(60)
  a <- s %o% c(1, 2, .5); b <- s %o% c(.8, 1.5, 2.5)
  k <- data.frame(participant = rep(c("P1", "P2"), each = 30),
                  trial = rep(1:3, each = 10, times = 2), cycle = rep(1:10, 6))
  p <- sample(60)
  original <- list(emg = a, motion = b[p, ])
  z <- alignObservationBlocks(original, list(emg = k, motion = k[p, ]), names(k))
  expect_identical(z$blocks, list(emg = a, motion = b))
  expect_identical(z$observations, k)
  expect_equal(z$row_index$motion, match(seq_len(60), p))
  expect_length(z$excluded_rows$motion, 0)
  expect_identical(original$motion, b[p, ])
  reference <- jive(list(emg = a, motion = b), rank_joint = 1, rank_individual = 0)
  expect_equal(jive(z$blocks, rank_joint = 1, rank_individual = 0), reference)
  expect_equal(multipleFactorAnalysis(z$blocks), multipleFactorAnalysis(list(emg = a, motion = b)))
})

test_that("unequal sets require an explicit policy and retain exclusions", {
  a <- matrix(1:6, 3); b <- matrix(7:12, 3)
  keys <- list(a = data.frame(id = c("x", "y", "z")),
               b = data.frame(id = c("z", "y", "w")))
  expect_error(alignObservationBlocks(list(a = a, b = b), keys, "id"), "sets differ")
  z <- alignObservationBlocks(list(a = a, b = b), keys, "id", "intersection")
  expect_identical(z$observations$id, c("y", "z"))
  expect_identical(z$row_index, list(a = 2:3, b = c(2L, 1L)))
  expect_identical(z$excluded_rows, list(a = 1L, b = 3L))
  keys$b$id <- rep("w", 3)
  expect_error(alignObservationBlocks(list(a = a, b = b), keys, "id"), "Duplicate")
  keys$b$id <- c("u", "v", "w")
  expect_error(alignObservationBlocks(list(a = a, b = b), keys, "id", "intersection"), "No shared")
})

test_that("key types and malformed inputs cannot silently align", {
  blocks <- list(a = matrix(1:4, 2), b = matrix(5:8, 2))
  keys <- list(a = data.frame(id = 1:2), b = data.frame(id = 1:2))
  expect_error(alignObservationBlocks(blocks, rev(keys), "id"), "same names")
  keys$b$id <- as.character(keys$b$id)
  expect_error(alignObservationBlocks(blocks, keys, "id"), "types must match")
  keys$b$id <- c(1L, NA_integer_)
  expect_error(alignObservationBlocks(blocks, keys, "id"), "complete")
  keys$b$id <- factor(1:2)
  expect_error(alignObservationBlocks(blocks, keys, "id"), "plain")
  keys$b$id <- 1:2; blocks$b[1, 1] <- Inf
  expect_error(alignObservationBlocks(blocks, keys, "id"), "finite numeric")
})

test_that("key strings cannot collide through delimiter concatenation", {
  k <- data.frame(participant = c("a:b", "a"), cycle = c("c", "b:c"))
  a <- matrix(1:4, 2)
  out <- alignObservationBlocks(list(a = a, b = a[2:1, ]),
        list(a = k, b = k[2:1, ]), names(k))
  expect_identical(out$blocks$a, out$blocks$b)
})
