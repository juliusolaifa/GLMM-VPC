# Calibration of the exact VPC-equality LRT. Uses helpers in difference_lrt.R.
# The distance reference is conservative for the vertex Gaussian experiment;
# no global GLMM uniformity or endpoint-maximality theorem is assumed.

#' Conservative Gaussian distance calibration for equal VPCs
#' @param statistic Nonnegative observed LRT.
#' @param beta Fixed intercept and slope under the null.
#' @param xa,xb Distinct covariate values.
#' @param family GLMM family supported by glmm_vpc().
#' @param family_params Named family-specific parameters.
#' @param information Efficient 3-by-3 covariance information (Schur complement).
#' @param alpha Significance levels.
#' @return List with conservative upper-bound p-value, critical values, and geometry.
calibrate_vpc_difference_distance <- function(
    statistic, beta, xa, xb, family = "poisson",
    family_params = numeric(0), information, alpha = c(.10, .05, .01)) {
  if (!is.numeric(statistic) || length(statistic) != 1L ||
      !is.finite(statistic) || statistic < 0)
    stop("statistic must be finite and nonnegative.", call. = FALSE)
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 0 | alpha >= 1))
    stop("alpha must contain probabilities strictly between zero and one.", call. = FALSE)
  contrast <- vpc_difference_vertex_contrast(beta, xa, xb, family, family_params)
  B <- cone_whitening_matrix(information)
  rays <- B %*% vpc_difference_vertex_rays(contrast)
  cosine <- sum(rays[, 1] * rays[, 2]) /
    sqrt(sum(rays[, 1]^2) * sum(rays[, 2]^2))
  angle <- acos(max(-1, min(1, cosine)))
  if (!is.finite(angle) || angle <= 0 || angle >= pi)
    stop("Degenerate wedge opening angle.", call. = FALSE)
  weights <- c(df1 = angle/(2*pi), df2 = .5, df3 = (pi-angle)/(2*pi))
  tailfun <- function(t) {
    sum(weights * stats::pchisq(t, df = 1:3, lower.tail = FALSE))
  }
  critical <- vapply(alpha, function(a) {
    upper <- 1
    while (tailfun(upper) > a) {
      upper <- 2 * upper
      if (!is.finite(upper)) stop("Failed to bracket critical value.")
    }
    stats::uniroot(function(t) tailfun(t) - a, c(0, upper), tol = 1e-10)$root
  }, numeric(1))
  names(critical) <- paste0("alpha_", format(alpha, trim = TRUE))
  p <- tailfun(statistic)
  list(statistic = statistic, p_value_upper_bound = p,
       critical_values = critical, reject = statistic > critical,
       angle = angle, weights = weights, contrast = contrast,
       method = "vertex_gaussian_conservative_distance",
       scope = "Gaussian vertex-local null displacements only")
}

#' Calibrate an equal-VPC LRT in a specified covariance regime
#'
#' Regime selection is explicit: a numerical rank threshold does not establish
#' uniform validity across vertex/rank-one/interior transitions.
#' @param regime One of vertex, rank_one, positive_definite.
#' @param information Vertex: efficient 3-by-3 information; otherwise full information.
#' @param gradient Full-parameter VPC-difference gradient for nonvertex regimes.
#' @param Gamma Nonzero rank-one covariance matrix for rank_one.
#' @param covariance_indices Covariance coordinate locations in full information.
#' @param nsim,seed Gaussian simulation settings for rank_one.
#' @param vertex_method Use the corrected H_Q bound by default; optional older distance bound.
#' @param ... Vertex arguments beta, xa, xb, family, family_params; or additional
#'   arguments accepted by the corresponding calibration routine.
calibrate_difference_lrt <- function(
    statistic, regime = c("vertex", "rank_one", "positive_definite"),
    information = NULL, gradient = NULL, Gamma = NULL,
    covariance_indices = 1:3, alpha = c(.10, .05, .01),
    nsim = 10000L, seed = NULL,
    vertex_method = c("improved_H", "distance"), ...) {
  regime <- match.arg(regime)
  vertex_method <- match.arg(vertex_method)
  if (!is.numeric(statistic) || length(statistic) != 1L ||
      !is.finite(statistic) || statistic < 0)
    stop("statistic must be a finite nonnegative scalar.")
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 0 | alpha >= 1))
    stop("alpha must contain probabilities strictly between zero and one.")
  if (regime == "vertex") {
    if (is.null(information)) stop("Vertex information is required.")
    fn <- if (vertex_method == "improved_H")
      calibrate_vpc_difference_improved else calibrate_vpc_difference_distance
    return(fn(statistic = statistic, information = information, alpha = alpha, ...))
  }
  if (is.null(information) || is.null(gradient))
    stop("Full information and gradient are required outside the vertex.")
  result <- calibrate_vpc_difference_regular_lrt(
    statistic = statistic, information = information, gradient = gradient,
    regime = if (regime == "rank_one") "rank_one" else "interior",
    Gamma = Gamma, covariance_indices = covariance_indices,
    nsim = nsim, seed = seed)
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 0 | alpha >= 1)) stop("Invalid alpha.")
  if (regime == "positive_definite") {
    result$critical_values <- stats::qchisq(1-alpha, df = 1)
    result$reject <- statistic > result$critical_values
  } else {
    result$critical_values <- vapply(alpha, function(a)
      as.numeric(stats::quantile(result$simulated, probs = 1-a,
                                 names = FALSE, type = 1)), numeric(1))
    result$reject <- statistic > result$critical_values
  }
  result$regime <- regime
  result$scope <- "Fixed-null asymptotic calibration; not globally uniform"
  result
}


# ---- Improved vertex Gaussian envelope: Theorem 6 / Appendix I9 ----
# IMPORTANT: The conservative Gaussian bound is NOT a globally uniform GLMM
# theorem, and it is NOT an alternative to maximizing the observed null likelihood.

.vpc_difference_canonical_Q <- function(contrast, information) {
  if (!is.matrix(information) || !identical(dim(information), c(3L,3L)) ||
      any(!is.finite(information)) ||
      max(abs(information-t(information))) > 1e-8 ||
      min(eigen(information, symmetric=TRUE, only.values=TRUE)$values) <= 0)
    stop("information must be a symmetric positive definite 3x3 matrix.")
  rays <- vpc_difference_vertex_rays(contrast)
  # The ray matrices are uu' with gamma ordering (11,12,22).
  uv <- lapply(1:2, function(i) {
    r <- rays[,i]
    u1 <- sqrt(max(0,r[1]))
    u2 <- if (u1 > 1e-12) r[2]/u1 else
      sign(r[2] + 1e-15)*sqrt(max(0,r[3]))
    c(u1,u2)
  })
  cross <- c(2*uv[[1]][1]*uv[[2]][1],
             uv[[1]][1]*uv[[2]][2] + uv[[2]][1]*uv[[1]][2],
             2*uv[[1]][2]*uv[[2]][2])
  L <- cbind(rays[,1],rays[,2],cross)
  if (abs(det(L)) < 1e-12*max(1,max(abs(L))^3))
    stop("Canonical covariance transformation is singular.")
  Q <- crossprod(L, information %*% L)
  list(Q=Q, L=L, rays=rays)
}

.vpc_orthant_distance <- function(x, A) {
  # Squared distance from x to A^(1/2) R_+^2 in Euclidean coordinates.
  S <- chol(A)
  # chol(A) is a valid square root for distance purposes.
  candidates <- list(c(0,0))
  for (j in 1:2) {
    g <- S[,j]
    candidates[[length(candidates)+1L]] <- g*max(0,sum(x*g)/sum(g*g))
  }
  ab <- tryCatch(solve(S,x), error=function(e) c(-1,-1))
  if (all(ab>=0)) candidates[[length(candidates)+1L]] <- x
  sqrt(min(vapply(candidates,function(y) sum((x-y)^2),numeric(1))))
}

.vpc_difference_envelope_geometry <- function(Q) {
  # Q coordinates are (a,b,z), as in manuscript equation (123).
  A <- Q[1:2,1:2,drop=FALSE]
  h <- Q[1:2,3]
  D <- as.numeric(Q[3,3]-crossprod(h,solve(A,h)))
  if (!is.finite(D) || D <= 0) stop("Nonpositive Schur complement.")
  k2 <- as.numeric(crossprod(h,solve(A,h))/D)
  k <- sqrt(max(0,k2))
  d <- as.vector(chol(A) %*% solve(A,h))/sqrt(D)
  lp <- .vpc_orthant_distance(d,A)
  lm <- .vpc_orthant_distance(-d,A)
  cosangle <- A[1,2]/sqrt(A[1,1]*A[2,2])
  angle <- acos(max(-1,min(1,cosangle)))
  list(kappa=k, M=1+k2, lambda_plus=lp,lambda_minus=lm,
       angle=angle, weights=c(angle/(2*pi),.5,(pi-angle)/(2*pi)))
}

#' Improved conservative Gaussian-local equality calibration (Theorem 6)
#'
#' Uses H_Q in closed form and optionally J_Q by deterministic numerical
#' integration. For numerical robustness, integration errors must be checked
#' before using the minimum as a conservative bound.
#' @export
calibrate_vpc_difference_improved <- function(
    statistic, beta, xa, xb, family="poisson",
    family_params=numeric(0), information, alpha=c(.05),
    use_J=FALSE, rel.tol=2e-5) {
  if (!is.logical(use_J) || length(use_J) != 1L || is.na(use_J))
    stop("use_J must be TRUE or FALSE.")
  if (!is.numeric(rel.tol) || length(rel.tol) != 1L ||
      !is.finite(rel.tol) || rel.tol <= 0 || rel.tol >= 1)
    stop("rel.tol must be a finite number in (0, 1).")
  if (!is.numeric(statistic) || length(statistic)!=1L ||
      !is.finite(statistic) || statistic<0)
    stop("Invalid statistic.")
  if (!is.numeric(alpha) || !length(alpha) ||
      any(!is.finite(alpha)) || any(alpha<=0 | alpha>=1))
    stop("Invalid alpha.")
  contrast <- vpc_difference_vertex_contrast(beta,xa,xb,family,family_params)
  can <- .vpc_difference_canonical_Q(contrast,information)
  g <- .vpc_difference_envelope_geometry(can$Q)
  p1 <- function(t) stats::pchisq(t,1,lower.tail=FALSE)
  H <- function(t) {
    if (t<=0) return(1)
    k <- g$kappa; M <- g$M
    .5*p1(t)+.5*p1(t/M)+
      (pi-g$angle)*k/(2*pi*sqrt(M))*exp(-t/(2*M))*
      (2*stats::pnorm(k*sqrt(t/M))-1)
  }
  # Numerical J_Q from equation (125) of vpc.pdf. Integration is
  # approximate, not a rigorous interval-certified upper bound.
  # The exact zero-coupling identity avoids a flat radial objective.
  J <- function(t) {
    if (t <= 0) return(1)
    L <- g$kappa
    if (L == 0) return(p1(t))
    M <- 1 + L * L
    halfnorm <- function(v) sqrt(2 / pi) * exp(-v * v / 2)
    radial_tail <- function(v, lam, df) {
      if (v * v >= t) return(1)
      delta <- t - v * v
      r0 <- sqrt(delta)
      if (r0 <= lam * v) {
        root <- r0
      } else {
        # Rationalized positive root of
        # L^2*r^2 + 2*lam*v*r = M*(t-v^2) + lam^2*v^2.
        c0 <- M * delta + (lam * v)^2
        root <- c0 / (sqrt((lam * v)^2 + L^2 * c0) + lam * v)
      }
      stats::pchisq(root * root, df = df, lower.tail = FALSE)
    }
    integrate_tail <- function(lam, df) {
      f <- function(v) vapply(v, function(vi)
        radial_tail(vi, lam, df) * halfnorm(vi), numeric(1))
      # Split at sqrt(t), where the radial tail becomes identically one.
      cutoff <- sqrt(t)
      integral <- stats::integrate(f, 0, cutoff, rel.tol = rel.tol,
                                   subdivisions = 500L,
                                   stop.on.error = TRUE)
      if (!is.finite(integral$value) || !is.finite(integral$abs.error))
        stop("J_Q quadrature failed.")
      # Upper numerical allowance: conservative relative to the returned
      # quadrature estimate, but NOT a mathematically certified error bound.
      estimate <- integral$value + stats::pchisq(t, 1, lower.tail = FALSE)
      min(1, estimate + 10 * integral$abs.error)
    }
    per_sign <- vapply(c(g$lambda_plus, g$lambda_minus), function(lam) {
      # The chi_0 component has R=0 and B=V^2.
      g$weights[1] * p1(t) +
        g$weights[2] * integrate_tail(lam, 1) +
        g$weights[3] * integrate_tail(lam, 2)
    }, numeric(1))
    min(1, max(0, mean(per_sign)))
  }
  tailfun <- function(t) {
    h <- H(t)
    if (!use_J) return(h)
    min(h,J(t))
  }
  critical <- vapply(alpha,function(a) {
    hi <- 1
    while (tailfun(hi)>a) {
      hi <- 2*hi
      if (!is.finite(hi)) stop("Failed to bracket critical value.")
    }
    stats::uniroot(function(t) tailfun(t)-a,c(0,hi),
                   tol=1e-6)$root
  },numeric(1))
  names(critical) <- paste0("alpha_",alpha)
  list(statistic=statistic,p_value_upper_bound=tailfun(statistic),
       critical_values=critical,reject=statistic>critical,
       geometry=g,canonical_information=can$Q,
       method=if(use_J) "min(J_Q_numeric,H_Q)" else "H_Q",
       scope="Gaussian vertex-local conservative; not globally uniform GLMM",
       numerical_integration=if(use_J) "approximate quadrature; not certified for strict size control" else "closed form")
}
