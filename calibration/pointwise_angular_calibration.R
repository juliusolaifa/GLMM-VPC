# Deterministic angular calibration for the pointwise GLMM VPC LRT.
# Requires the existing helpers: project_whitened_psd(),
# whitened_pointwise_ray(), and .cm_information().
# No regime argument, random draws, or RNG state changes.

.pw_gauss_legendre <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n != floor(n) || n < 2L) stop('n must be an integer >= 2')
  n <- as.integer(n)
  k <- seq_len(n - 1L)
  a <- k / sqrt(4 * k * k - 1)
  A <- matrix(0, n, n)
  A[cbind(k, k + 1L)] <- a
  A[cbind(k + 1L, k)] <- a
  ee <- eigen(A, symmetric = TRUE)
  ord <- order(ee$values)
  list(nodes = ee$values[ord], weights = 2 * ee$vectors[1L, ord]^2)
}

# Precompute a deterministic spherical quadrature rule and cone projections.
# The rule is Gauss-Legendre in cos(polar angle) and periodic trapezoid
# in azimuth; weights sum to one for uniform U on S^2.
pointwise_angular_grid <- function(x0, information, n_polar = 20L,
                                   n_azimuth = 48L,
                                   projection_intervals = 64L,
                                   projection_tol = 1e-10) {
  J <- .cm_information(information)
  if (!is.numeric(x0) || length(x0) != 1L || !is.finite(x0))
    stop('x0 must be a finite numeric scalar')
  .check_count <- function(n, minimum, name) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
        n != floor(n) || n < minimum || n > .Machine$integer.max)
      stop(sprintf('%s must be an integer >= %d', name, minimum))
    as.integer(n)
  }
  n_polar <- .check_count(n_polar, 2L, 'n_polar')
  n_azimuth <- .check_count(n_azimuth, 4L, 'n_azimuth')
  projection_intervals <- .check_count(projection_intervals, 4L, 'projection_intervals')
  if (!is.numeric(projection_tol) || length(projection_tol) != 1L ||
      !is.finite(projection_tol) || projection_tol <= 0)
    stop('projection_tol must be a positive finite scalar')
  gl <- .pw_gauss_legendre(n_polar)
  az <- 2 * pi * (seq_len(n_azimuth) - 0.5) / n_azimuth
  e <- as.numeric(whitened_pointwise_ray(x0, J))
  e <- e / sqrt(sum(e * e))
  N <- n_polar * n_azimuth
  d <- numeric(N)
  weight <- numeric(N)
  idx <- 0L
  for (i in seq_len(n_polar)) {
    z <- gl$nodes[i]
    r <- sqrt(max(0, 1 - z * z))
    for (j in seq_len(n_azimuth)) {
      idx <- idx + 1L
      u <- c(r * cos(az[j]), r * sin(az[j]), z)
      proj <- project_whitened_psd(u, J,
                  intervals = projection_intervals,
                  tol = projection_tol)$projection
      dc <- sum(proj^2) - max(0, sum(u * e))^2
      if (dc < -1e-7) stop('Negative projection difference: inspect geometry')
      d[idx] <- max(0, dc)
      weight[idx] <- gl$weights[i] / (2 * n_azimuth)
    }
  }
  structure(list(d = d, weight = weight, x0 = x0,
                 information = J, n_polar = n_polar,
                 n_azimuth = n_azimuth), class = 'pointwise_angular_grid')
}

# Vertex survival and density at t > 0, using R^2 ~ chisq(3).
pointwise_vertex_survival <- function(t, grid) {
  stopifnot(inherits(grid, 'pointwise_angular_grid'))
  if (!is.numeric(t) || any(!is.finite(t)) || any(t < 0))
    stop('t must be finite and nonnegative')
  vapply(t, function(q) {
    if (q == 0) return(sum(grid$weight[grid$d > 0]))
    ok <- grid$d > 0
    sum(grid$weight[ok] * stats::pchisq(q / grid$d[ok], 3,
                                        lower.tail = FALSE))
  }, numeric(1))
}

pointwise_vertex_density <- function(t, grid) {
  stopifnot(inherits(grid, 'pointwise_angular_grid'))
  if (!is.numeric(t) || any(!is.finite(t)) || any(t <= 0))
    stop('t must be finite and positive')
  ok <- grid$d > 0
  vapply(t, function(q)
    sum(grid$weight[ok] * stats::dchisq(q / grid$d[ok], 3) /
        grid$d[ok]), numeric(1))
}

pointwise_ray_survival <- function(t) {
  0.5 * stats::pchisq(t, 1, lower.tail = FALSE) +
    0.5 * stats::pchisq(t, 2, lower.tail = FALSE)
}
pointwise_ray_density <- function(t) {
  0.5 * stats::dchisq(t, 1) + 0.5 * stats::dchisq(t, 2)
}

# Single envelope calibration for any pointwise-null fit.
calibrate_pointwise_angular <- function(statistic, x0, information,
                                        alpha = 0.05, grid = NULL,
                                        n_polar = 20L, n_azimuth = 48L,
                                        projection_intervals = 64L) {
  if (length(statistic) != 1L || !is.finite(statistic) || statistic < 0)
    stop('statistic must be finite and nonnegative')
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop('alpha must be in (0,1)')
  J <- .cm_information(information)
  if (!is.numeric(x0) || length(x0) != 1L || !is.finite(x0))
    stop('x0 must be a finite numeric scalar')
  if (is.null(grid)) {
    grid <- pointwise_angular_grid(
      x0, J, n_polar, n_azimuth, projection_intervals)
  } else {
    if (!inherits(grid, 'pointwise_angular_grid') ||
        !isTRUE(all.equal(grid$x0, x0, tolerance = 1e-12)) ||
        !isTRUE(all.equal(grid$information, J, tolerance = 1e-10)))
      stop('grid must match the supplied x0 and information')
    if (!is.numeric(grid$d) || !is.numeric(grid$weight) ||
        length(grid$d) != length(grid$weight) ||
        any(!is.finite(grid$d)) || any(grid$d < 0) ||
        any(!is.finite(grid$weight)) || any(grid$weight < 0) ||
        !isTRUE(all.equal(sum(grid$weight), 1, tolerance = 1e-8)))
      stop('invalid angular grid weights or projection values')
  }
  pv <- pointwise_vertex_survival(statistic, grid)
  pr <- pointwise_ray_survival(statistic)
  quantile_tail <- function(fun) {
    # For a distribution with an atom at zero, its upper-tail critical
    # value can be zero at sufficiently large alpha.
    if (fun(0) <= alpha) return(0)
    hi <- 8
    for (k in seq_len(100L)) {
      if (fun(hi) <= alpha) break
      hi <- hi * 2
    }
    if (!is.finite(hi) || fun(hi) > alpha)
      stop('unable to bracket angular critical value')
    stats::uniroot(function(q) fun(q) - alpha,
                   c(0, hi), tol = 1e-8)$root
  }
  cv <- quantile_tail(function(q) pointwise_vertex_survival(q, grid))
  cr <- quantile_tail(pointwise_ray_survival)
  list(statistic = statistic, x0 = x0, alpha = alpha,
       p_value = max(pv, pr), critical_value = max(cv, cr),
       reject = statistic > max(cv, cr),
       vertex_p_value = pv, ray_p_value = pr,
       vertex_critical = cv, ray_critical = cr,
       method = 'deterministic_angular_envelope', grid = grid)
}

# Independent refinement check. Compare outputs from two angular grids.
pointwise_angular_refinement <- function(x0, information,
                                         thresholds = c(1, 3, 5, 8),
                                         coarse = c(12L, 24L),
                                         fine = c(24L, 48L),
                                         projection_intervals = 64L) {
  if (length(coarse) != 2L || length(fine) != 2L)
    stop('coarse and fine must each contain two grid dimensions')
  g1 <- pointwise_angular_grid(x0, information, coarse[1], coarse[2],
                               projection_intervals)
  g2 <- pointwise_angular_grid(x0, information, fine[1], fine[2],
                               projection_intervals)
  p1 <- pointwise_vertex_survival(thresholds, g1)
  p2 <- pointwise_vertex_survival(thresholds, g2)
  data.frame(threshold = thresholds, coarse_tail = p1,
             fine_tail = p2, absolute_difference = abs(p1 - p2))
}
