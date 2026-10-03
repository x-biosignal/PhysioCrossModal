# Nonparametric bivariate Granger causality via Wilson spectral factorization.

.gc_scalar_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop(sprintf("'%s' must be one integer of at least %d", name, minimum),
         call. = FALSE)
  }
  as.integer(x)
}

.gc_default_segment_length <- function(n, overlap) {
  candidates <- 2^(seq.int(4L, floor(log2(n))))
  usable <- candidates + pmax(1L, floor(candidates * (1 - overlap))) <= n
  candidates <- candidates[usable]
  if (!length(candidates)) {
    stop("Signals are too short for two spectral realizations", call. = FALSE)
  }
  as.integer(min(512L, max(candidates)))
}

.gc_segment_starts <- function(n, segment_length, overlap) {
  step <- max(1L, as.integer(floor(segment_length * (1 - overlap))))
  seq.int(1L, n - segment_length + 1L, by = step)
}

.gc_outer_spectrum <- function(x, y, taper, n_fft, sr) {
  padded_x <- complex(length.out = n_fft)
  padded_y <- complex(length.out = n_fft)
  n <- length(x)
  padded_x[seq_len(n)] <- x * taper
  padded_y[seq_len(n)] <- y * taper
  fourier <- cbind(stats::fft(padded_x), stats::fft(padded_y))
  normalization <- sr * sum(taper^2)
  out <- array(0 + 0i, dim = c(2L, 2L, n_fft))
  for (frequency in seq_len(n_fft)) {
    z <- fourier[frequency, ]
    out[, , frequency] <- z %o% Conj(z) / normalization
  }
  out
}

.estimate_gc_csd <- function(
    x, y, sr, spectral_estimator, n_fft, segment_length, overlap,
    time_bandwidth, n_tapers) {
  n <- length(x)
  if (!is.numeric(overlap) || length(overlap) != 1L ||
      !is.finite(overlap) || overlap < 0 || overlap >= 1) {
    stop("'overlap' must be one finite number in [0, 1)", call. = FALSE)
  }
  if (is.null(segment_length)) {
    segment_length <- .gc_default_segment_length(n, overlap)
  } else {
    segment_length <- .gc_scalar_integer(
      segment_length, "segment_length", minimum = 16L
    )
    if (segment_length > n) {
      stop("'segment_length' cannot exceed signal length", call. = FALSE)
    }
  }
  if (is.null(n_fft)) {
    n_fft <- as.integer(2^ceiling(log2(segment_length)))
  } else {
    n_fft <- .gc_scalar_integer(n_fft, "n_fft", minimum = 16L)
  }
  if (n_fft %% 2L != 0L || n_fft < segment_length) {
    stop("'n_fft' must be even and at least 'segment_length'", call. = FALSE)
  }

  starts <- .gc_segment_starts(n, segment_length, overlap)
  if (identical(spectral_estimator, "welch") && length(starts) < 2L) {
    stop("Welch estimation requires at least two complete segments",
         call. = FALSE)
  }

  if (identical(spectral_estimator, "multitaper")) {
    if (!is.numeric(time_bandwidth) || length(time_bandwidth) != 1L ||
        !is.finite(time_bandwidth) || time_bandwidth <= 1 ||
        time_bandwidth >= segment_length / 2) {
      stop(
        "'time_bandwidth' must be finite, greater than 1, and below half ",
        "the segment length",
        call. = FALSE
      )
    }
    if (is.null(n_tapers)) {
      n_tapers <- as.integer(floor(2 * time_bandwidth) - 1L)
    } else {
      n_tapers <- .gc_scalar_integer(n_tapers, "n_tapers", minimum = 2L)
    }
    if (n_tapers < 2L || n_tapers > segment_length) {
      stop("'n_tapers' must be between 2 and 'segment_length'",
           call. = FALSE)
    }
    tapers <- .dpss_tapers(segment_length, time_bandwidth, n_tapers)
  } else {
    tapers <- matrix(.window_hanning(segment_length), ncol = 1L)
    n_tapers <- 1L
  }

  spectrum <- array(0 + 0i, dim = c(2L, 2L, n_fft))
  n_realizations <- 0L
  for (start in starts) {
    index <- start:(start + segment_length - 1L)
    segment_x <- x[index] - mean(x[index])
    segment_y <- y[index] - mean(y[index])
    for (taper_index in seq_len(ncol(tapers))) {
      spectrum <- spectrum + .gc_outer_spectrum(
        segment_x, segment_y, tapers[, taper_index], n_fft, sr
      )
      n_realizations <- n_realizations + 1L
    }
  }
  if (n_realizations < 2L) {
    stop("Spectral estimation requires at least two realizations",
         call. = FALSE)
  }
  spectrum <- spectrum / n_realizations

  n_one_sided <- n_fft %/% 2L + 1L
  spectrum <- spectrum[, , seq_len(n_one_sided), drop = FALSE]
  if (n_one_sided > 2L) {
    spectrum[, , 2:(n_one_sided - 1L)] <-
      2 * spectrum[, , 2:(n_one_sided - 1L), drop = FALSE]
  }
  for (frequency in seq_len(n_one_sided)) {
    current <- spectrum[, , frequency]
    spectrum[, , frequency] <- (current + Conj(t(current))) / 2
  }

  list(
    spectrum = spectrum,
    frequencies = (0:(n_one_sided - 1L)) * sr / n_fft,
    n_fft = n_fft,
    segment_length = segment_length,
    n_realizations = n_realizations,
    n_tapers = n_tapers
  )
}

.regularize_gc_spectrum <- function(spectrum, eigen_floor) {
  if (!is.numeric(eigen_floor) || length(eigen_floor) != 1L ||
      !is.finite(eigen_floor) || eigen_floor <= 0 || eigen_floor >= 1) {
    stop("'eigen_floor' must be one finite number in (0, 1)",
         call. = FALSE)
  }

  maximum_adjustment <- 0
  for (frequency in seq_len(dim(spectrum)[3L])) {
    original <- spectrum[, , frequency]
    original <- (original + Conj(t(original))) / 2
    scale <- mean(Re(diag(original)))
    if (!is.finite(scale) || scale <= 0) {
      stop("Cross-spectral matrix has non-positive power", call. = FALSE)
    }
    decomposition <- eigen(original, symmetric = TRUE)
    values <- Re(decomposition$values)
    material_negative <- max(100 * eigen_floor, 1e-8) * scale
    if (min(values) < -material_negative) {
      stop(
        sprintf(
          "Cross-spectral matrix is materially indefinite at bin %d",
          frequency
        ),
        call. = FALSE
      )
    }
    floor_value <- eigen_floor * scale
    adjusted_values <- pmax(values, floor_value)
    adjusted <- decomposition$vectors %*% diag(adjusted_values, 2L) %*%
      Conj(t(decomposition$vectors))
    adjusted <- (adjusted + Conj(t(adjusted))) / 2
    relative <- sqrt(sum(Mod(adjusted - original)^2)) /
      max(sqrt(sum(Mod(original)^2)), .Machine$double.eps)
    maximum_adjustment <- max(maximum_adjustment, relative)
    spectrum[, , frequency] <- adjusted
  }

  list(spectrum = spectrum, maximum_adjustment = maximum_adjustment)
}

.gc_one_to_two_spectrum <- function(spectrum) {
  n_one_sided <- dim(spectrum)[3L]
  n_fft <- 2L * (n_one_sided - 1L)
  output <- array(0 + 0i, dim = c(2L, 2L, n_fft))
  output[, , 1L] <- 2 * spectrum[, , 1L]
  if (n_one_sided > 2L) {
    for (frequency in 2:(n_one_sided - 1L)) {
      output[, , frequency] <- spectrum[, , frequency]
      negative <- n_fft - frequency + 2L
      output[, , negative] <- t(spectrum[, , frequency])
    }
  }
  output[, , n_one_sided] <- 2 * spectrum[, , n_one_sided]
  output
}

.fft_matrix_array <- function(x, inverse = FALSE) {
  output <- array(0 + 0i, dim = dim(x))
  n_fft <- dim(x)[3L]
  for (row in seq_len(dim(x)[1L])) {
    for (column in seq_len(dim(x)[2L])) {
      transformed <- stats::fft(x[row, column, ], inverse = inverse)
      if (inverse) {
        transformed <- transformed / n_fft
      }
      output[row, column, ] <- transformed
    }
  }
  output
}

.wilson_plus <- function(x) {
  n_fft <- dim(x)[3L]
  if (length(dim(x)) != 3L || dim(x)[1L] != dim(x)[2L] ||
      n_fft < 4L || n_fft %% 2L != 0L || any(!is.finite(x))) {
    stop("Wilson plus input must be a finite square matrix spectrum on an ",
         "even grid",
         call. = FALSE)
  }

  lag <- .fft_matrix_array(x, inverse = TRUE)
  positive <- array(0 + 0i, dim = dim(x))
  zero_lag <- lag[, , 1L]
  zero_lag[lower.tri(zero_lag)] <- 0
  diag(zero_lag) <- diag(zero_lag) / 2
  positive[, , 1L] <- zero_lag
  last_positive <- n_fft %/% 2L
  positive[, , 2:last_positive] <- lag[, , 2:last_positive, drop = FALSE]
  .fft_matrix_array(positive)
}

.wilson_sfactor <- function(
    spectral_matrix, tolerance = 1e-8, max_iterations = 1000L) {
  if (length(dim(spectral_matrix)) != 3L ||
      dim(spectral_matrix)[1L] != dim(spectral_matrix)[2L] ||
      dim(spectral_matrix)[1L] < 2L ||
      dim(spectral_matrix)[3L] < 4L ||
      dim(spectral_matrix)[3L] %% 2L != 0L ||
      any(!is.finite(spectral_matrix))) {
    stop("'spectral_matrix' must be a finite square matrix spectrum on an ",
         "even two-sided grid",
         call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1) {
    stop("'tolerance' must be one finite number in (0, 1)",
         call. = FALSE)
  }
  max_iterations <- .gc_scalar_integer(
    max_iterations, "max_iterations", minimum = 1L
  )

  n_channel <- dim(spectral_matrix)[1L]
  n_fft <- dim(spectral_matrix)[3L]
  identity <- diag(n_channel)
  structural_tolerance <- min(
    1e-6, max(100 * tolerance, 100 * .Machine$double.eps)
  )
  for (frequency in seq_len(n_fft)) {
    current <- spectral_matrix[, , frequency]
    hermitian_error <- max(Mod(current - Conj(t(current))))
    if (hermitian_error > structural_tolerance * max(1, max(Mod(current)))) {
      stop("Wilson input is not Hermitian", call. = FALSE)
    }
  }

  covariance <- Re(.fft_matrix_array(spectral_matrix, inverse = TRUE)[, , 1L])
  covariance <- (covariance + t(covariance)) / 2
  initial <- tryCatch(
    chol(covariance),
    error = function(error) {
      stop("Wilson initialization failed: lag-zero covariance is not ",
           "positive definite",
           call. = FALSE)
    }
  )
  psi <- array(0 + 0i, dim = dim(spectral_matrix))
  for (frequency in seq_len(n_fft)) {
    psi[, , frequency] <- initial
  }

  converged <- FALSE
  relative_update <- Inf
  iterations <- 0L
  for (iteration in seq_len(max_iterations)) {
    g <- array(0 + 0i, dim = dim(spectral_matrix))
    for (frequency in seq_len(n_fft)) {
      inverse_psi <- tryCatch(
        solve(psi[, , frequency]),
        error = function(error) {
          stop("Wilson iteration produced a singular spectral factor",
               call. = FALSE)
        }
      )
      current <- inverse_psi %*% spectral_matrix[, , frequency] %*%
        Conj(t(inverse_psi)) + identity
      g[, , frequency] <- (current + Conj(t(current))) / 2
    }
    plus <- .wilson_plus(g)
    updated <- array(0 + 0i, dim = dim(psi))
    for (frequency in seq_len(n_fft)) {
      updated[, , frequency] <- psi[, , frequency] %*%
        plus[, , frequency]
    }
    relative_update <- sqrt(sum(Mod(updated - psi)^2)) /
      max(sqrt(sum(Mod(psi)^2)), .Machine$double.eps)
    psi <- updated
    iterations <- as.integer(iteration)
    if (is.finite(relative_update) && relative_update < tolerance) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    stop(
      sprintf(
        paste0(
          "Wilson factorization did not converge after %d iterations ",
          "(relative update %.3g)"
        ),
        max_iterations, relative_update
      ),
      call. = FALSE
    )
  }

  reconstructed <- array(0 + 0i, dim = dim(spectral_matrix))
  for (frequency in seq_len(n_fft)) {
    reconstructed[, , frequency] <- psi[, , frequency] %*%
      Conj(t(psi[, , frequency]))
  }
  reconstruction_error <- sqrt(sum(Mod(reconstructed - spectral_matrix)^2)) /
    max(sqrt(sum(Mod(spectral_matrix)^2)), .Machine$double.eps)
  reconstruction_limit <- min(0.05, max(0.02, 100 * tolerance))
  if (!is.finite(reconstruction_error) ||
      reconstruction_error > reconstruction_limit) {
    stop(
      sprintf(
        "Wilson spectral reconstruction error %.3g exceeds %.3g",
        reconstruction_error, reconstruction_limit
      ),
      call. = FALSE
    )
  }

  lag_factor <- .fft_matrix_array(psi, inverse = TRUE)
  a0 <- lag_factor[, , 1L]
  imaginary_scale <- max(abs(Im(a0)))
  real_scale <- max(1, max(abs(Re(a0))))
  if (imaginary_scale > structural_tolerance * real_scale) {
    stop("Wilson zero-lag factor is not real within tolerance",
         call. = FALSE)
  }
  a0 <- Re(a0)
  sigma <- a0 %*% t(a0)
  sigma <- (sigma + t(sigma)) / 2
  if (any(!is.finite(sigma)) || any(diag(sigma) <= 0)) {
    stop("Wilson innovation covariance is not positive on the diagonal",
         call. = FALSE)
  }
  inverse_a0 <- tryCatch(
    solve(a0),
    error = function(error) {
      stop("Wilson zero-lag factor is singular", call. = FALSE)
    }
  )
  transfer <- array(0 + 0i, dim = dim(psi))
  for (frequency in seq_len(n_fft)) {
    transfer[, , frequency] <- psi[, , frequency] %*% inverse_a0
  }

  list(
    psi = psi,
    H = transfer,
    Sigma = sigma,
    reconstructed = reconstructed,
    converged = converged,
    iterations = iterations,
    relative_update = relative_update,
    reconstruction_error = reconstruction_error,
    reconstruction_limit = reconstruction_limit
  )
}

.gc_trapezoid_mean <- function(values, frequencies, range) {
  selected <- which(frequencies >= range[1L] & frequencies <= range[2L])
  if (length(selected) < 2L) {
    stop("'freq_range' must contain at least two spectral bins",
         call. = FALSE)
  }
  selected_frequencies <- frequencies[selected]
  selected_values <- values[selected]
  span <- selected_frequencies[length(selected_frequencies)] -
    selected_frequencies[1L]
  if (!is.finite(span) || span <= 0) {
    stop("'freq_range' has no positive frequency span", call. = FALSE)
  }
  left <- seq_len(length(selected_values) - 1L)
  right <- left + 1L
  sum(diff(selected_frequencies) *
        (selected_values[left] + selected_values[right]) / 2) / span
}

.gc_nonparametric <- function(
    x, y, sr, freq_range, spectral_estimator, n_fft, segment_length, overlap,
    time_bandwidth, n_tapers, factor_tolerance, factor_max_iterations,
    eigen_floor) {
  spectral_estimator <- match.arg(
    spectral_estimator, c("multitaper", "welch")
  )
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) ||
      length(x) < 32L || any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Nonparametric Granger requires equal-length finite numeric signals ",
         "with at least 32 samples",
         call. = FALSE)
  }
  if (!is.numeric(sr) || length(sr) != 1L || !is.finite(sr) || sr <= 0) {
    stop("'sr' must be one positive finite sampling rate", call. = FALSE)
  }
  if (!is.numeric(factor_tolerance) || length(factor_tolerance) != 1L ||
      !is.finite(factor_tolerance) || factor_tolerance <= 0 ||
      factor_tolerance >= 1) {
    stop("'factor_tolerance' must be one finite number in (0, 1)",
         call. = FALSE)
  }
  factor_max_iterations <- .gc_scalar_integer(
    factor_max_iterations, "factor_max_iterations", minimum = 1L
  )
  if (!is.null(freq_range) &&
      (!is.numeric(freq_range) || length(freq_range) != 2L ||
       any(!is.finite(freq_range)) || freq_range[1L] < 0 ||
       freq_range[1L] >= freq_range[2L] || freq_range[2L] > sr / 2)) {
    stop("'freq_range' must be two ordered finite values inside [0, sr / 2]",
         call. = FALSE)
  }

  x <- as.numeric(x) - mean(x)
  y <- as.numeric(y) - mean(y)
  scale_x <- sqrt(mean(x^2))
  scale_y <- sqrt(mean(y^2))
  if (!is.finite(scale_x) || !is.finite(scale_y) ||
      scale_x <= .Machine$double.eps || scale_y <= .Machine$double.eps) {
    stop("Nonparametric Granger requires non-constant signals",
         call. = FALSE)
  }
  normalized_x <- x / scale_x
  normalized_y <- y / scale_y
  relation <- sum(normalized_x * normalized_y) /
    sqrt(sum(normalized_x^2) * sum(normalized_y^2))
  residual <- normalized_y - relation * normalized_x
  if (abs(relation) > 1 - 1e-12 &&
      sqrt(mean(residual^2)) < 1e-8) {
    stop("Signals are rank-deficient or numerically identical; spectral ",
         "factorization is undefined",
         call. = FALSE)
  }

  estimate <- .estimate_gc_csd(
    x, y, sr, spectral_estimator, n_fft, segment_length, overlap,
    time_bandwidth, n_tapers
  )
  regularized <- .regularize_gc_spectrum(estimate$spectrum, eigen_floor)
  two_sided <- .gc_one_to_two_spectrum(regularized$spectrum)
  factor <- .wilson_sfactor(
    two_sided,
    tolerance = factor_tolerance,
    max_iterations = factor_max_iterations
  )

  n_one_sided <- length(estimate$frequencies)
  reconstructed <- factor$reconstructed[
    , , seq_len(n_one_sided), drop = FALSE
  ]
  transfer <- factor$H[, , seq_len(n_one_sided), drop = FALSE]
  sigma <- factor$Sigma
  conditional_y <- sigma[2L, 2L] - sigma[1L, 2L]^2 / sigma[1L, 1L]
  conditional_x <- sigma[1L, 1L] - sigma[1L, 2L]^2 / sigma[2L, 2L]
  if (!is.finite(conditional_x) || !is.finite(conditional_y) ||
      conditional_x <= 0 || conditional_y <= 0) {
    stop("Conditional innovation variances are not positive",
         call. = FALSE)
  }

  power_x <- Re(reconstructed[1L, 1L, ])
  power_y <- Re(reconstructed[2L, 2L, ])
  denominator_yx <- power_x - conditional_y * Mod(transfer[1L, 2L, ])^2
  denominator_xy <- power_y - conditional_x * Mod(transfer[2L, 1L, ])^2
  denominator_tolerance <- min(
    0.05, max(100 * factor_tolerance, 1e-8)
  )
  if (any(denominator_yx <= -denominator_tolerance * power_x) ||
      any(denominator_xy <= -denominator_tolerance * power_y)) {
    stop("Geweke decomposition produced a materially invalid denominator",
         call. = FALSE)
  }
  denominator_yx <- pmax(denominator_yx, .Machine$double.eps * power_x)
  denominator_xy <- pmax(denominator_xy, .Machine$double.eps * power_y)
  ratio_yx <- power_x / denominator_yx
  ratio_xy <- power_y / denominator_xy
  if (any(ratio_yx < 1 - denominator_tolerance) ||
      any(ratio_xy < 1 - denominator_tolerance)) {
    stop("Geweke decomposition produced materially negative causality",
         call. = FALSE)
  }
  gc_yx_spectrum <- pmax(0, log(pmax(1, ratio_yx)))
  gc_xy_spectrum <- pmax(0, log(pmax(1, ratio_xy)))
  if (any(!is.finite(gc_xy_spectrum)) ||
      any(!is.finite(gc_yx_spectrum))) {
    stop("Geweke causality spectrum is not finite", call. = FALSE)
  }

  aggregation_range <- if (is.null(freq_range)) {
    c(0, sr / 2)
  } else {
    freq_range
  }
  gc_xy <- .gc_trapezoid_mean(
    gc_xy_spectrum, estimate$frequencies, aggregation_range
  )
  gc_yx <- .gc_trapezoid_mean(
    gc_yx_spectrum, estimate$frequencies, aggregation_range
  )

  list(
    gc_xy = gc_xy,
    gc_yx = gc_yx,
    net_gc = gc_xy - gc_yx,
    order = NA_integer_,
    method = "nonparametric",
    frequencies = estimate$frequencies,
    gc_xy_spectrum = gc_xy_spectrum,
    gc_yx_spectrum = gc_yx_spectrum,
    spectral_estimator = spectral_estimator,
    converged = factor$converged,
    iterations = factor$iterations,
    factorization_error = factor$reconstruction_error,
    factorization_limit = factor$reconstruction_limit,
    regularization_adjustment = regularized$maximum_adjustment,
    innovation_covariance = sigma,
    n_realizations = estimate$n_realizations,
    n_fft = estimate$n_fft,
    segment_length = estimate$segment_length,
    n_tapers = estimate$n_tapers,
    freq_range = aggregation_range
  )
}
