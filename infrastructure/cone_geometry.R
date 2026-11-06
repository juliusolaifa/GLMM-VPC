# Shared internal validation and covariance-coordinate helpers.
.cg_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop(sprintf("%s must be a finite numeric scalar.", name), call. = FALSE)
  as.numeric(x)
}

.cg_vector3 <- function(x, name) {
  if (!is.numeric(x) || length(x) != 3L || any(!is.finite(x)))
    stop(sprintf("%s must be a finite numeric vector of length three.", name),
         call. = FALSE)
  as.numeric(x)
}

# Covariance coordinates: (gamma11, gamma12, gamma22).
gamma_to_matrix <- function(gamma) {
  g <- .cg_vector3(gamma, "gamma")
  matrix(c(g[1L], g[2L], g[2L], g[3L]), nrow = 2L)
}


matrix_to_gamma <- function(Gamma) {
  if (!is.matrix(Gamma) || !is.numeric(Gamma) ||
      !identical(dim(Gamma), c(2L, 2L)) || any(!is.finite(Gamma)) ||
      !isTRUE(all.equal(Gamma, t(Gamma), tolerance = 1e-12))) {
    stop("Gamma must be a finite symmetric 2-by-2 matrix.", call. = FALSE)
  }
  c(gamma11 = Gamma[1L, 1L], gamma12 = (Gamma[1L, 2L] + Gamma[2L, 1L])/2,
    gamma22 = Gamma[2L, 2L])
}

#' Linear covariance-to-Lorentz coordinate transformation
#' @return The 3-by-3 matrix T satisfying u = T %*% gamma.
#' @export
lorentz_coordinate_matrix <- function() {
  s <- 1 / sqrt(2)
  rbind(c(s, 0, s), c(s, 0, -s), c(0, sqrt(2), 0))
}

#' Convert covariance coordinates to Lorentz coordinates
#' @param gamma Numeric covariance vector of length three.
#' @return Numeric Lorentz coordinate vector of length three.
#' @export
gamma_to_lorentz <- function(gamma) {
  drop(lorentz_coordinate_matrix() %*% .cg_vector3(gamma, "gamma"))
}

#' Convert Lorentz coordinates to covariance coordinates
#' @param u Numeric Lorentz vector of length three.
#' @return Named covariance vector.
#' @export
lorentz_to_gamma <- function(u) {
  g <- drop(solve(lorentz_coordinate_matrix(), .cg_vector3(u, "u")))
  stats::setNames(g, c("gamma11", "gamma12", "gamma22"))
}

#' Project onto the standard three-dimensional Lorentz cone
#' @param u Numeric vector (t, y1, y2).
#' @return Euclidean projection onto t >= sqrt(y1^2+y2^2).
#' @export
project_lorentz <- function(u) {
  u <- .cg_vector3(u, "u")
  t <- u[1L]
  y <- u[2:3]
  r <- sqrt(sum(y * y))
  if (r <= t) return(u)
  if (r <= -t) return(c(0, 0, 0))
  a <- (t + r) / 2
  c(a, a * y / r)
}

#' Project a symmetric matrix onto the positive-semidefinite cone
#' @param Gamma Finite symmetric 2-by-2 matrix.
#' @return Frobenius-norm projection onto the PSD cone.
#' @export
project_psd2 <- function(Gamma) {
  g <- matrix_to_gamma(Gamma)
  M <- gamma_to_matrix(g)
  E <- eigen(M, symmetric = TRUE)
  P <- E$vectors %*% diag(pmax(E$values, 0), nrow = 2L) %*% t(E$vectors)
  (P + t(P)) / 2
}

#' Project onto a fixed nonnegative ray
#' @param x Finite numeric vector.
#' @param direction Finite nonzero numeric vector of the same length.
#' @return Euclidean projection onto {t * direction : t >= 0}.
#' @export
project_ray <- function(x, direction) {
  if (!is.numeric(x) || !is.numeric(direction) ||
      !length(x) || length(x) != length(direction) ||
      any(!is.finite(x)) || any(!is.finite(direction)) ||
      sum(direction^2) <= 0) {
    stop("x and direction must be finite equal-length vectors; direction must be nonzero.",
         call. = FALSE)
  }
  pmax(0, sum(x * direction) / sum(direction^2)) * direction
}

#' Null ray for zero random-effects variance at x0
#' @param x0 Finite covariate value.
#' @param coordinates One of 'gamma', 'lorentz', or 'matrix'.
#' @return A generator of the PSD null ray Gamma %*% c(1,x0) = 0.
#' @export
pointwise_null_ray <- function(x0, coordinates = c("gamma", "lorentz", "matrix")) {
  x0 <- .cg_scalar(x0, "x0")
  coordinates <- match.arg(coordinates)
  # Orthogonal to z(x0) = (1, x0): q = (-x0, 1).
  q <- c(-x0, 1)
  q <- q / sqrt(sum(q^2))
  M <- tcrossprod(q)
  if (coordinates == "matrix") return(M)
  g <- matrix_to_gamma(M)
  if (coordinates == "gamma") return(g)
  gamma_to_lorentz(g)
}

#' Compute Euclidean projection-distance difference for nested cones
#' @param observation Numeric vector.
#' @param projection_null Its projection onto the null cone.
#' @param projection_full Its projection onto the full cone.
#' @return Squared-distance difference (null minus full).
#' @export
projection_lrt_difference <- function(observation, projection_null, projection_full) {
  x <- observation
  n <- length(x)
  valid <- function(v) is.numeric(v) && length(v) == n && all(is.finite(v))
  if (n == 0L || !valid(x) || !valid(projection_null) || !valid(projection_full)) {
    stop("All inputs must be finite numeric vectors of the same nonzero length.",
         call. = FALSE)
  }
  sum((x - projection_null)^2) - sum((x - projection_full)^2)
}


# Efficient-information metric geometry for 2x2 covariance parameters.
# Coordinates: gamma = (gamma11, gamma12, gamma22).

.cm_information <- function(information) {
  if (is.list(information) && !is.null(information$I_metric))
    information <- information$I_metric
  if (!is.matrix(information) || !is.numeric(information) ||
      !identical(dim(information), c(3L, 3L)) ||
      any(!is.finite(information)) ||
      !isTRUE(all.equal(information, t(information), tolerance = 1e-12)))
    stop("information must be a finite symmetric 3-by-3 matrix or an object containing I_metric.", call. = FALSE)
  information <- (information + t(information))/2
  if (inherits(try(chol(information), silent = TRUE), "try-error"))
    stop("information must be positive definite.", call. = FALSE)
  information
}

#' Efficient-information whitening map
#'
#' Returns a matrix B such that ||B %*% gamma||^2 = gamma' I gamma.
#' The upper-triangular Cholesky factor is used, so crossproducts are exact.
#' @param information Positive-definite 3-by-3 efficient information matrix,
#'   or a list with an I_metric element.
#' @return A 3-by-3 matrix B with crossprod(B) = information.
#' @export
cone_whitening_matrix <- function(information) {
  chol(.cm_information(information))
}

#' Transform covariance coordinates into information-whitened coordinates
#' @param gamma Covariance vector (gamma11, gamma12, gamma22).
#' @param information Efficient information matrix.
#' @return Numeric vector of length three.
#' @export
whiten_gamma <- function(gamma, information) {
  drop(cone_whitening_matrix(information) %*% .cg_vector3(gamma, "gamma"))
}

#' Transform whitened coordinates back to covariance coordinates
#' @param y Numeric whitened vector of length three.
#' @param information Efficient information matrix.
#' @return Named covariance-coordinate vector.
#' @export
unwhiten_gamma <- function(y, information) {
  g <- drop(backsolve(cone_whitening_matrix(information), .cg_vector3(y, "y")))
  stats::setNames(g, c("gamma11", "gamma12", "gamma22"))
}

#' Pointwise null ray in information-whitened coordinates
#' @param x0 Covariate value at which the random-effects variance is zero.
#' @param information Efficient information matrix.
#' @param unit Normalize the returned direction to Euclidean length one.
#' @return A length-three ray generator.
#' @export
whitened_pointwise_ray <- function(x0, information, unit = TRUE) {
  if (!is.logical(unit) || length(unit) != 1L || is.na(unit))
    stop("unit must be TRUE or FALSE.", call. = FALSE)
  direction <- whiten_gamma(pointwise_null_ray(x0, "gamma"), information)
  if (unit) direction <- direction / sqrt(sum(direction^2))
  direction
}

#' Project a whitened observation onto a pointwise null ray
#' @param y A whitened length-three observation.
#' @param x0 Null covariate value.
#' @param information Efficient information matrix.
#' @return Euclidean projection in whitened coordinates.
#' @export
project_pointwise_ray <- function(y, x0, information) {
  project_ray(.cg_vector3(y, "y"), whitened_pointwise_ray(x0, information))
}

#' Information-metric projection onto a pointwise null ray
#' @param gamma A covariance-coordinate vector, not necessarily PSD.
#' @param x0 Null covariate value.
#' @param information Efficient information matrix.
#' @return Covariance-coordinate projection onto the null ray.
#' @export
project_pointwise_ray_metric <- function(gamma, x0, information) {
  I <- .cm_information(information)
  direction <- unname(pointwise_null_ray(x0, "gamma"))
  g <- .cg_vector3(gamma, "gamma")
  coefficient <- max(0, drop(crossprod(direction, I %*% g)) /
                       drop(crossprod(direction, I %*% direction)))
  stats::setNames(coefficient * direction,
                  c("gamma11", "gamma12", "gamma22"))
}


# Projection onto the information-whitened 2x2 PSD cone.
# Requires cone_geometry.R and cone_metric.R.

#' Project onto the information-transformed PSD cone
#'
#' Computes the Euclidean projection of y onto
#' {B gamma : gamma_to_matrix(gamma) is PSD}, where
#' crossprod(B) equals the supplied efficient information matrix.
#' The boundary is parameterized by rank-one PSD matrices, and a
#' one-dimensional search is used to find the global projection.
#'
#' @param y Finite numeric vector of length three, in whitened coordinates.
#' @param information Positive-definite 3-by-3 efficient information matrix,
#'   or a list containing I_metric.
#' @param intervals Number of angular intervals searched (at least 32).
#' @param tol Absolute tolerance for each bounded scalar optimization.
#' @return A list with projection, gamma, squared_distance, and region.
#' @export
project_whitened_psd <- function(y, information, intervals = 128L,
                                 tol = 1e-11) {
  y <- .cg_vector3(y, "y")
  I <- .cm_information(information)
  if (!is.numeric(intervals) || length(intervals) != 1L ||
      !is.finite(intervals) || intervals != as.integer(intervals) ||
      intervals < 32L) {
    stop("intervals must be an integer >= 32.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive finite number.", call. = FALSE)
  }
  B <- chol(I)
  g <- drop(backsolve(B, y))
  eig <- eigen(gamma_to_matrix(g), symmetric = TRUE, only.values = TRUE)$values
  # Interior and boundary points already in the transformed cone.
  if (min(eig) >= 0) {
    return(list(projection = y, gamma = stats::setNames(g,
                c("gamma11", "gamma12", "gamma22")),
                squared_distance = 0, region = "inside"))
  }

  # Every nonzero boundary point is t * q(theta) q(theta)' with t >= 0.
  # For fixed theta, eliminate t analytically via ray projection.
  # Search every angular interval, not only local maxima on a coarse grid.
  score <- function(theta) {
    a <- cos(theta)
    b <- sin(theta)
    d <- drop(B %*% c(a*a, a*b, b*b))
    h <- sum(y*d)
    if (h <= 0) return(0)
    h*h / sum(d*d)
  }
  grid <- seq(0, pi, length.out = as.integer(intervals) + 1L)
  values <- vapply(grid, score, numeric(1))
  best_index <- which.max(values)
  best_theta <- grid[best_index]
  best_score <- values[best_index]
  for (j in seq_len(length(grid) - 1L)) {
    candidate <- stats::optimize(function(t) -score(t),
                                 interval = grid[c(j, j + 1L)], tol = tol)
    candidate_score <- -candidate$objective
    if (candidate_score > best_score) {
      best_score <- candidate_score
      best_theta <- candidate$minimum
    }
  }
  a <- cos(best_theta)
  b <- sin(best_theta)
  d_gamma <- c(a*a, a*b, b*b)
  d <- drop(B %*% d_gamma)
  coefficient <- max(0, sum(y*d) / sum(d*d))
  projection <- coefficient * d
  gamma <- coefficient * d_gamma
  list(projection = projection,
       gamma = stats::setNames(gamma, c("gamma11", "gamma12", "gamma22")),
       squared_distance = sum((y - projection)^2),
       region = if (coefficient == 0) "vertex" else "boundary")
}

#' Project covariance coordinates onto the PSD cone in Fisher metric
#' @param gamma Finite covariance-coordinate vector of length three.
#' @param information Positive-definite efficient information matrix.
#' @param ... Passed to project_whitened_psd().
#' @return A list with projected gamma, whitened projection, distance, region.
#' @export
project_psd_metric <- function(gamma, information, ...) {
  y <- whiten_gamma(gamma, information)
  project_whitened_psd(y, information, ...)
}

#' Projection-difference LRT for a pointwise zero-variance null
#' @param y Whitened Gaussian observation.
#' @param x0 Covariate value of the pointwise null.
#' @param information Efficient information matrix.
#' @param ... Passed to project_whitened_psd().
#' @return Nonnegative squared-distance difference, subject to rounding.
#' @export
pointwise_projection_statistic <- function(y, x0, information, ...) {
  y <- .cg_vector3(y, "y")
  full <- project_whitened_psd(y, information, ...)$projection
  null <- project_pointwise_ray(y, x0, information)
  statistic <- projection_lrt_difference(y, null, full)
  if (statistic < -1e-8 * max(1, sum(y*y))) {
    stop("Projection nesting check failed; inspect numerical accuracy.",
         call. = FALSE)
  }
  max(0, statistic)
}
