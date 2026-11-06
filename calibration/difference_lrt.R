# Geometry for equality of random-effects variances at two covariate values.
# Coordinates are (gamma11, gamma12, gamma22).
# IMPORTANT: equality of v(x) is NOT generally equality of GLMM VPC rho(x).

#' Variance-difference contrast
#'
#' Returns c such that v(xa)-v(xb) = crossprod(c, gamma), where
#' v(x) = gamma11 + 2*x*gamma12 + x^2*gamma22.
#' @param xa,xb Distinct finite covariate values.
#' @return Named length-three numeric vector.
#' @export
variance_difference_contrast <- function(xa, xb) {
  xa <- .cg_scalar(xa, "xa")
  xb <- .cg_scalar(xb, "xb")
  if (xa == xb) stop("xa and xb must be distinct.", call. = FALSE)
  stats::setNames(c(0, 2 * (xa - xb), xa^2 - xb^2),
                  c("gamma11", "gamma12", "gamma22"))
}

#' Extreme rays of the equal-variance PSD wedge
#'
#' The null set {gamma: Gamma >= 0, v(xa)=v(xb)} is the cone generated
#' by the two rank-one PSD directions returned as columns.
#' @param xa,xb Distinct finite covariate values.
#' @return 3-by-2 matrix in covariance-vector coordinates.
#' @export
variance_equality_rays <- function(xa, xb) {
  variance_difference_contrast(xa, xb)
  m <- (xa + xb) / 2
  rays <- cbind(c(1, 0, 0), c(m*m, -m, 1))
  rownames(rays) <- c("gamma11", "gamma12", "gamma22")
  colnames(rays) <- c("constant_intercept", "balanced_predictor")
  rays
}

#' Project onto a two-generator closed convex wedge
#'
#' Minimizes ||y - D a||^2 over a >= 0. The exact candidate set
#' consists of the vertex, the two extreme rays, and the unconstrained
#' two-dimensional projection when its coefficients are nonnegative.
#' @param y Finite length-three vector.
#' @param directions Full-column-rank 3-by-2 matrix of wedge generators.
#' @return List with projection, coefficients and squared_distance.
#' @export
project_two_ray_wedge <- function(y, directions) {
  y <- .cg_vector3(y, "y")
  D <- directions
  if (!is.matrix(D) || !is.numeric(D) || !identical(dim(D), c(3L, 2L)) ||
      any(!is.finite(D)) || qr(D)$rank != 2L)
    stop("directions must be a finite, rank-two 3-by-2 matrix.", call. = FALSE)
  candidates <- list(c(0, 0))
  for (j in 1:2) {
    a <- c(0, 0)
    a[j] <- max(0, sum(y * D[, j]) / sum(D[, j]^2))
    candidates[[length(candidates) + 1L]] <- a
  }
  interior <- drop(solve(crossprod(D), crossprod(D, y)))
  if (all(interior >= 0)) candidates[[length(candidates) + 1L]] <- interior
  distances <- vapply(candidates, function(a) sum((y - drop(D %*% a))^2), numeric(1))
  best <- candidates[[which.min(distances)]]
  list(projection = drop(D %*% best), coefficients = best,
       squared_distance = min(distances))
}

#' Information-whitened equal-variance wedge
#' @param xa,xb Distinct finite covariate values.
#' @param information Positive-definite 3-by-3 efficient information matrix.
#' @return Matrix whose columns generate the whitened null wedge.
#' @export
whitened_variance_equality_wedge <- function(xa, xb, information) {
  B <- cone_whitening_matrix(information)
  B %*% variance_equality_rays(xa, xb)
}

#' Gaussian-limit equal-variance LRT at the covariance vertex
#'
#' Computes distance(y,C0)^2 - distance(y,C)^2, where C is the
#' information-whitened PSD cone and C0 is the equal-variance wedge.
#' This is a variance-equality statistic, NOT a general VPC-equality LRT.
#' @param y Whitened Gaussian observation, length three.
#' @param xa,xb Distinct covariate values.
#' @param information Positive-definite efficient information matrix.
#' @param ... Passed to project_whitened_psd().
#' @return Nonnegative numeric statistic (up to projection tolerance).
#' @export
variance_difference_lrt_limit_statistic <- function(y, xa, xb, information, ...) {
  y <- .cg_vector3(y, "y")
  D <- whitened_variance_equality_wedge(xa, xb, information)
  null_dist <- project_two_ray_wedge(y, D)$squared_distance
  full_dist <- project_whitened_psd(y, information, ...)$squared_distance
  stat <- null_dist - full_dist
  if (stat < -1e-8 * max(1, sum(y*y)))
    stop("Projection nesting check failed; inspect numerical accuracy.", call. = FALSE)
  max(0, stat)
}

#' Simulate the vertex-null equal-variance LRT limit
#'
#' Uses Gaussian projection-difference simulation; no chi-bar-square
#' assumption is imposed on the difference statistic.
#' @param xa,xb Distinct covariate values.
#' @param information Efficient information matrix.
#' @param nsim Positive integer number of simulations.
#' @param seed Optional nonnegative integer; caller RNG state is restored.
#' @param ... Passed to project_whitened_psd().
#' @return Numeric vector of limiting statistics.
#' @export
simulate_variance_difference_vertex_lrt <- function(xa, xb, information,
                                                     nsim = 10000L, seed = NULL, ...) {
  variance_difference_contrast(xa, xb)
  I <- .cm_information(information)
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != floor(nsim) || nsim > .Machine$integer.max)
    stop("nsim must be a positive integer.", call. = FALSE)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed != floor(seed) || seed > .Machine$integer.max)
      stop("seed must be a nonnegative integer.", call. = FALSE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  draws <- matrix(stats::rnorm(3L * as.integer(nsim)), ncol = 3L)
  vapply(seq_len(as.integer(nsim)), function(i)
    variance_difference_lrt_limit_statistic(draws[i, ], xa, xb, I, ...), numeric(1))
}

#' Monte Carlo calibration of a vertex equal-variance LRT
#'
#' This function applies ONLY when the scientific null is equivalent to
#' v(xa)=v(xb). For a general GLMM VPC difference rho(xa)=rho(xb),
#' the null also involves the fixed effects and observation variance;
#' a separate model-specific constrained likelihood is required.
#' @param statistic Nonnegative observed LRT statistic.
#' @param xa,xb Distinct covariate values.
#' @param information Efficient information matrix.
#' @param nsim,seed Monte Carlo controls.
#' @param ... Projection numerical controls.
#' @return List with p_value, statistic, method, simulated and mcse.
#' @export
calibrate_variance_difference_vertex_lrt <- function(statistic, xa, xb,
                                                      information, nsim = 10000L,
                                                      seed = NULL, ...) {
  statistic <- .cg_scalar(statistic, "statistic")
  if (statistic < 0) stop("statistic must be nonnegative.", call. = FALSE)
  simulated <- simulate_variance_difference_vertex_lrt(xa, xb, information,
                                                        nsim, seed, ...)
  p <- (1 + sum(simulated >= statistic)) / (length(simulated) + 1)
  list(statistic = statistic, p_value = p,
       method = "equal_variance_vertex_projection_monte_carlo",
       nsim = length(simulated), simulated = simulated,
       mcse = sqrt(p * (1 - p) / (length(simulated) + 1)))
}


# VPC-equality constraint for random-intercept/slope GLMMs.
# This module evaluates the exact scientific null; it does not fit a
# restricted GLMM or claim a universal difference-LRT calibration.

#' VPC difference at two covariate values
#' @param beta Fixed intercept and slope.
#' @param Gamma Symmetric PSD random-effect covariance matrix.
#' @param xa,xb Distinct finite covariate values.
#' @param family Response family supported by glmm_vpc().
#' @param family_params Family-specific parameters.
#' @return Scalar rho(xa) - rho(xb).
#' @export
vpc_difference <- function(beta, Gamma, xa, xb,
                           family = c("poisson", "nbinom", "cpg", "gaussian"),
                           family_params = numeric(0)) {
  family <- match.arg(family)
  if (!is.numeric(xa) || length(xa) != 1L || !is.finite(xa) ||
      !is.numeric(xb) || length(xb) != 1L || !is.finite(xb) || xa == xb) {
    stop("xa and xb must be distinct finite scalars.", call. = FALSE)
  }
  rho <- glmm_vpc(beta, Gamma, c(xa, xb), family, family_params)
  unname(rho[1L] - rho[2L])
}

#' Vertex derivative of VPC with respect to random-effect variance
#'
#' At Gamma=0, rho(x) = k(x) * v(x) + o(||Gamma||), where
#' v(x) = (1,x) Gamma (1,x)^T.
#' @inheritParams vpc_difference
#' @param x Finite numeric covariate values.
#' @return Positive numeric coefficients k(x).
#' @export
vpc_vertex_scale <- function(beta, x,
                             family = c("poisson", "nbinom", "cpg", "gaussian"),
                             family_params = numeric(0)) {
  family <- match.arg(family)
  if (!is.numeric(beta) || length(beta) != 2L || any(!is.finite(beta)) ||
      !is.numeric(x) || any(!is.finite(x))) {
    stop("beta must have length two and x must be finite numeric.", call. = FALSE)
  }
  # Reuse glmm_vpc's family validation at the vertex.
  invisible(glmm_vpc(beta, matrix(0, 2, 2), x, family, family_params))
  eta <- beta[1L] + beta[2L] * x
  if (family == "gaussian") return(rep(1 / unname(family_params["sigma2"]), length(x)))
  if (family == "poisson") return(exp(eta))
  if (family == "nbinom") return(1 / (exp(-eta) + 1 / unname(family_params["theta"])))
  exp((2 - unname(family_params["power"])) * eta) / unname(family_params["phi"])
}

#' Linearized VPC-difference null at the covariance vertex
#'
#' For gamma=(gamma11,gamma12,gamma22), the derivative of
#' rho(xa)-rho(xb) at Gamma=0 is crossprod(c, gamma).
#' Unlike variance_difference_contrast(), c generally has a nonzero
#' gamma11 coefficient. Nuisance derivatives vanish at Gamma=0.
#' @inheritParams vpc_difference
#' @return Named length-three gradient in covariance-vector coordinates.
#' @export
vpc_difference_vertex_contrast <- function(beta, xa, xb,
                                           family = c("poisson", "nbinom", "cpg", "gaussian"),
                                           family_params = numeric(0)) {
  family <- match.arg(family)
  if (!is.numeric(xa) || length(xa) != 1L || !is.finite(xa) ||
      !is.numeric(xb) || length(xb) != 1L || !is.finite(xb) || xa == xb)
    stop("xa and xb must be distinct finite scalars.", call. = FALSE)
  k <- vpc_vertex_scale(beta, c(xa, xb), family, family_params)
  q <- function(x) c(1, 2*x, x*x)
  stats::setNames(k[1L]*q(xa) - k[2L]*q(xb),
                  c("gamma11", "gamma12", "gamma22"))
}

#' Check the exact VPC-equality null
#' @inheritParams vpc_difference
#' @param tolerance Nonnegative numerical tolerance.
#' @return List with difference, satisfies_null and tolerance.
#' @export
check_vpc_difference_null <- function(beta, Gamma, xa, xb,
                                      family = c("poisson", "nbinom", "cpg", "gaussian"),
                                      family_params = numeric(0), tolerance = 1e-8) {
  family <- match.arg(family)
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0)
    stop("tolerance must be a finite nonnegative scalar.", call. = FALSE)
  d <- vpc_difference(beta, Gamma, xa, xb, family, family_params)
  list(difference = d, satisfies_null = abs(d) <= tolerance, tolerance = tolerance)
}


# Local Gaussian-limit calibration for equality of GLMM VPCs at Gamma = 0.
# Requires cone_geometry.R, cone_metric.R, cone_projection.R,
# difference_lrt.R, and vpc_difference_constraint.R.

#' Extreme rays of the linearized VPC-equality cone at the covariance vertex
#'
#' Given c = gradient_gamma{rho(xa)-rho(xb)} at Gamma=0, the local
#' null is {Gamma >= 0: c' vech(Gamma)=0}. For distinct xa and xb
#' and positive vertex VPC scales, this is a two-ray convex wedge.
#' @param contrast Numeric length-three vector ordered (gamma11,gamma12,gamma22).
#' @param tolerance Positive numerical threshold for eigenvalue checks.
#' @return 3-by-2 matrix of rank-one PSD generators.
#' @export
vpc_difference_vertex_rays <- function(contrast, tolerance = 1e-10) {
  if (!is.numeric(contrast) || length(contrast) != 3L ||
      any(!is.finite(contrast)))
    stop("contrast must be a finite numeric vector of length three.", call. = FALSE)
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0)
    stop("tolerance must be a positive finite scalar.", call. = FALSE)
  # For rank-one Gamma = (u,v)(u,v)', c'gamma = (u,v)' A (u,v).
  A <- matrix(c(contrast[1L], contrast[2L] / 2,
                contrast[2L] / 2, contrast[3L]), 2L, 2L)
  ee <- eigen(A, symmetric = TRUE)
  scale <- max(abs(ee$values))
  if (scale == 0 || min(abs(ee$values)) <= tolerance * scale ||
      prod(ee$values) >= 0)
    stop("contrast does not define a numerically nondegenerate two-ray wedge.",
         call. = FALSE)
  ip <- which.max(ee$values)
  ineg <- which.min(ee$values)
  vp <- ee$vectors[, ip]
  vn <- ee$vectors[, ineg]
  a <- sqrt(-ee$values[ineg])
  b <- sqrt(ee$values[ip])
  to_gamma <- function(t) c(t[1L]^2, t[1L] * t[2L], t[2L]^2)
  rays <- cbind(to_gamma(a * vp + b * vn),
                to_gamma(a * vp - b * vn))
  rownames(rays) <- c("gamma11", "gamma12", "gamma22")
  colnames(rays) <- c("ray1", "ray2")
  rays
}

#' Gaussian-limit VPC-difference statistic at the covariance vertex
#' @param y Whitened Gaussian score vector of length three.
#' @param information Positive definite 3-by-3 efficient information matrix.
#' @param contrast Linearized VPC-difference contrast at Gamma=0.
#' @param ... Arguments forwarded to project_whitened_psd().
#' @return Nonnegative projection-difference statistic.
#' @export
vpc_difference_vertex_lrt_limit_statistic <- function(y, information, contrast, ...) {
  y <- .cg_vector3(y, "y")
  B <- cone_whitening_matrix(information)
  D <- B %*% vpc_difference_vertex_rays(contrast)
  d0 <- project_two_ray_wedge(y, D)$squared_distance
  d1 <- project_whitened_psd(y, information, ...)$squared_distance
  stat <- d0 - d1
  if (stat < -1e-8 * max(1, sum(y^2)))
    stop("Projection nesting check failed; inspect numerical accuracy.", call. = FALSE)
  max(0, stat)
}

#' Monte Carlo calibration of the VPC-difference LRT at Gamma=0
#'
#' Uses the local cone generated by the VPC difference derivative, not
#' the equal-variance cone. The result is a vertex-specific asymptotic
#' calibration, not a finite-sample or uniform boundary guarantee.
#' @param statistic Observed nonnegative LRT statistic.
#' @param beta Fixed intercept and slope at the null vertex.
#' @param xa,xb Distinct covariate values.
#' @param family GLMM family.
#' @param family_params Family-specific parameters.
#' @param information Efficient covariance information at the null vertex.
#' @param nsim Positive integer number of Gaussian draws.
#' @param seed Optional RNG seed; caller RNG state is restored.
#' @param ... Passed to the full-cone projection.
#' @return List containing p_value, mcse, simulated, contrast and method.
#' @export
calibrate_vpc_difference_vertex_lrt <- function(
    statistic, beta, xa, xb,
    family = c("poisson", "nbinom", "cpg", "gaussian"),
    family_params = numeric(0), information,
    nsim = 10000L, seed = NULL, ...) {
  if (!is.numeric(statistic) || length(statistic) != 1L ||
      !is.finite(statistic) || statistic < 0)
    stop("statistic must be a finite nonnegative scalar.", call. = FALSE)
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != floor(nsim) || nsim > .Machine$integer.max)
    stop("nsim must be a positive integer.", call. = FALSE)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed != floor(seed) || seed > .Machine$integer.max)
      stop("seed must be a nonnegative integer.", call. = FALSE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  family <- match.arg(family)
  contrast <- vpc_difference_vertex_contrast(beta, xa, xb, family, family_params)
  # Fail early if geometry or information is degenerate.
  vpc_difference_vertex_rays(contrast)
  B <- cone_whitening_matrix(information)
  D <- B %*% vpc_difference_vertex_rays(contrast)
  draws <- matrix(stats::rnorm(3L * as.integer(nsim)), ncol = 3L)
  simulated <- vapply(seq_len(as.integer(nsim)), function(i) {
    y <- draws[i, ]
    d0 <- project_two_ray_wedge(y, D)$squared_distance
    d1 <- project_whitened_psd(y, information, ...)$squared_distance
    stat <- d0 - d1
    if (stat < -1e-8 * max(1, sum(y^2)))
      stop("Projection nesting check failed; inspect numerical accuracy.", call. = FALSE)
    max(0, stat)
  }, numeric(1))
  p <- (1 + sum(simulated >= statistic)) / (length(simulated) + 1)
  list(statistic = statistic, p_value = p,
       mcse = sqrt(p * (1 - p) / (length(simulated) + 1)),
       nsim = length(simulated), simulated = simulated,
       contrast = contrast, method = "vpc_difference_vertex_projection_monte_carlo")
}


# Gaussian tangent-cone calibration away from the covariance vertex.
# Coordinates: theta = (gamma11, gamma12, gamma22, nuisance...),
# or use covariance_indices to locate the covariance coordinates.

.vd_regular_setup <- function(information, gradient, regime,
                              Gamma = NULL, covariance_indices = 1:3,
                              tolerance = 1e-9) {
  if (!is.matrix(information) || !is.numeric(information) ||
      nrow(information) != ncol(information) || any(!is.finite(information)) ||
      max(abs(information - t(information))) > 1e-10)
    stop("information must be a finite symmetric matrix.", call. = FALSE)
  d <- nrow(information)
  if (!is.numeric(gradient) || length(gradient) != d || any(!is.finite(gradient)))
    stop("gradient must be a finite vector matching information.", call. = FALSE)
  if (!is.numeric(covariance_indices) || length(covariance_indices) != 3L ||
      anyNA(covariance_indices) || any(covariance_indices != as.integer(covariance_indices)) ||
      any(covariance_indices < 1L | covariance_indices > d) ||
      anyDuplicated(covariance_indices))
    stop("covariance_indices must identify three distinct coordinates.", call. = FALSE)
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0)
    stop("tolerance must be positive.", call. = FALSE)
  B <- tryCatch(chol(information), error = function(e) NULL)
  if (is.null(B)) stop("information must be positive definite.", call. = FALSE)
  g <- drop(forwardsolve(t(B), gradient)) # B^{-T} gradient
  ng <- sqrt(sum(g*g))
  if (ng <= tolerance) stop("VPC-difference gradient is degenerate.", call. = FALSE)
  g <- g/ng
  if (regime == "interior") return(list(g = g, n = NULL, regime = regime))
  if (is.null(Gamma) || !is.matrix(Gamma) || !identical(dim(Gamma), c(2L,2L)) ||
      !is.numeric(Gamma) || any(!is.finite(Gamma)) ||
      max(abs(Gamma-t(Gamma))) > tolerance)
    stop("Gamma must be a symmetric finite 2-by-2 matrix.", call. = FALSE)
  ev <- eigen((Gamma+t(Gamma))/2, symmetric = TRUE)
  scale <- max(1, max(abs(ev$values)))
  if (ev$values[1L] <= tolerance*scale || abs(ev$values[2L]) > tolerance*scale)
    stop("Gamma must be nonzero rank one PSD within tolerance.", call. = FALSE)
  u <- ev$vectors[, 2L]
  a <- numeric(d)
  a[covariance_indices] <- c(u[1L]^2, 2*u[1L]*u[2L], u[2L]^2)
  n <- drop(forwardsolve(t(B), a))
  n <- n/sqrt(sum(n*n))
  n_plane <- n - sum(n*g)*g
  if (sqrt(sum(n_plane*n_plane)) <= tolerance)
    stop("Equality gradient is parallel to the PSD boundary normal; nontransverse case requires separate theory.", call. = FALSE)
  list(g = g, n = n, regime = regime)
}

#' Tangent-cone Gaussian LRT for VPC equality away from Gamma=0
#'
#' The gradient must be for the EXACT VPC difference with respect to the
#' FULL parameter vector, including fixed effects and family parameters.
#' The information matrix must use the same parameter coordinates and include
#' nuisance parameters; do not substitute covariance-only information unless
#' the nuisance gradients are genuinely zero or nuisance is known.
#' @param z Standard Gaussian score in whitened coordinates.
#' @param information Positive definite information matrix for full parameter vector.
#' @param gradient Gradient of rho(xa)-rho(xb) in matching coordinates.
#' @param regime Either "interior" or "rank_one".
#' @param Gamma Nonzero rank-one 2-by-2 PSD matrix when regime="rank_one".
#' @param covariance_indices Indices of gamma11, gamma12, gamma22.
#' @param tolerance Numerical tolerance.
#' @return Nonnegative projection-difference statistic.
#' @export
vpc_difference_regular_limit_statistic <- function(
    z, information, gradient, regime = c("interior", "rank_one"),
    Gamma = NULL, covariance_indices = 1:3, tolerance = 1e-9) {
  regime <- match.arg(regime)
  s <- .vd_regular_setup(information, gradient, regime, Gamma,
                         covariance_indices, tolerance)
  if (!is.numeric(z) || length(z) != length(s$g) || any(!is.finite(z)))
    stop("z must be a finite whitened vector matching information.", call. = FALSE)
  # Full tangent: R^d in interior, {v: n'v >= 0} at rank-one boundary.
  full <- z
  if (!is.null(s$n)) full <- z + max(0, -sum(s$n*z))*s$n
  # Equality tangent: g'v=0, intersected with the halfspace at rank one.
  nul <- z - sum(s$g*z)*s$g
  if (!is.null(s$n)) {
    np <- s$n - sum(s$n*s$g)*s$g
    nul <- nul + max(0, -sum(s$n*nul))/sum(np*np)*np
  }
  stat <- sum((z-nul)^2) - sum((z-full)^2)
  if (stat < -tolerance*max(1, sum(z*z)))
    stop("Negative nested-projection statistic; inspect geometry.", call. = FALSE)
  max(0, stat)
}

#' Calibrate VPC-difference LRT away from the covariance vertex
#'
#' Interior: asymptotic chi-square(1). Rank one: Monte Carlo simulation of
#' the joint halfspace/equality tangent-cone projection difference.
#' This is a pointwise asymptotic calibration, not a uniform boundary result.
#' @param statistic Observed nonnegative LRT statistic.
#' @param nsim Number of Gaussian simulations for rank-one calibration.
#' @param seed Optional RNG seed; caller RNG state is restored.
#' @inheritParams vpc_difference_regular_limit_statistic
#' @return List containing p_value, method, statistic, and simulation diagnostics.
#' @export
calibrate_vpc_difference_regular_lrt <- function(
    statistic, information, gradient, regime = c("interior", "rank_one"),
    Gamma = NULL, covariance_indices = 1:3, nsim = 10000L,
    seed = NULL, tolerance = 1e-9) {
  regime <- match.arg(regime)
  if (!is.numeric(statistic) || length(statistic) != 1L ||
      !is.finite(statistic) || statistic < 0)
    stop("statistic must be a nonnegative finite scalar.", call. = FALSE)
  s <- .vd_regular_setup(information, gradient, regime, Gamma,
                         covariance_indices, tolerance)
  if (regime == "interior") {
    return(list(statistic = statistic, p_value = stats::pchisq(statistic, 1L,
                lower.tail = FALSE), method = "interior_chisq1", nsim = 0L,
                mcse = NA_real_, simulated = NULL))
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != floor(nsim) || nsim > .Machine$integer.max)
    stop("nsim must be a positive integer.", call. = FALSE)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed != floor(seed) || seed > .Machine$integer.max)
      stop("seed must be a nonnegative integer.", call. = FALSE)
    had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had) old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had) assign(".Random.seed", old, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  draws <- matrix(stats::rnorm(as.integer(nsim)*length(s$g)),
                  ncol = length(s$g))
  simulated <- vapply(seq_len(nrow(draws)), function(i)
    vpc_difference_regular_limit_statistic(draws[i, ], information, gradient,
      regime = "rank_one", Gamma = Gamma,
      covariance_indices = covariance_indices, tolerance = tolerance), numeric(1))
  p <- (1 + sum(simulated >= statistic))/(length(simulated)+1)
  list(statistic = statistic, p_value = p,
       mcse = sqrt(p*(1-p)/(length(simulated)+1)), nsim = length(simulated),
       simulated = simulated, method = "rank_one_tangent_projection_monte_carlo")
}

# Observed restricted fitting and vpc_difference_lrt() are implemented only
# in glmmtmb_difference_adapter.R. Do not reintroduce a second optimizer here.
# This file owns difference geometry and Gaussian-limit reference laws.
