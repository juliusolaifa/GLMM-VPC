# Interior delta-method inference for response-scale GLMM VPCs.
# Parameter order: beta0, beta1, family-specific psi, gamma11, gamma12, gamma22.
# The information matrix must be supplied in exactly this order.
#
# This file intentionally does not estimate nonvertex Fisher information:
# fim_glmm.R currently targets the covariance vertex for NB2/CPG.

.vpc_ci_params <- function(fit) {
  if (!inherits(fit, "glmmVPC_fit"))
    stop("`fit` must be a glmmVPC_fit object.", call. = FALSE)
  if (length(fit$beta) != 2L || any(!is.finite(fit$beta)) ||
      !is.matrix(fit$Gamma) || !identical(dim(fit$Gamma), c(2L, 2L)))
    stop("Invalid beta or Gamma in fit.", call. = FALSE)
  family <- match.arg(fit$family, c("gaussian", "poisson", "nbinom", "cpg"))
  psi <- fit$family_params
  expected <- switch(family, gaussian = "sigma2", poisson = character(),
                     nbinom = "theta", cpg = c("phi", "power"))
  if (!is.numeric(psi) || !setequal(names(psi), expected) ||
      length(psi) != length(expected) || any(!is.finite(psi)))
    stop("Invalid family_params names or values.", call. = FALSE)
  psi <- psi[expected]
  # Delegate family-domain and PSD checks to the package's VPC evaluator.
  invisible(glmm_vpc(fit$beta, fit$Gamma, 0, family, psi))
  list(beta = unname(fit$beta), Gamma = fit$Gamma,
       family = family, psi = psi,
       names = c("beta0", "beta1", expected,
                 "gamma11", "gamma12", "gamma22"))
}

#' Analytical full-parameter gradient of the response-scale VPC
#'
#' The ordering is beta0, beta1, family parameters, gamma11, gamma12,
#' gamma22. At v=0, derivatives are obtained directly from B/(B+W),
#' not from log(B). This is an ambient gradient where the extension exists.
#' @param fit A glmmVPC_fit.
#' @param x One finite covariate value.
#' @return Named numeric gradient.
vpc_gradient <- function(fit, x) {
  a <- .vpc_ci_params(fit)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop("`x` must be one finite number.", call. = FALSE)
  eta <- sum(a$beta * c(1, x))
  v <- as.numeric(c(1, x) %*% a$Gamma %*% c(1, x))
  if (v < -1e-8) stop("Negative random-effects variance.", call. = FALSE)
  v <- max(v, 0)
  rho <- as.numeric(glmm_vpc(a$beta, a$Gamma, x, a$family, a$psi))
  dr_eta <- 0
  dr_psi <- setNames(numeric(length(a$psi)), names(a$psi))
  if (v == 0) {
    # B_v(0)/W(0); B_eta(0)=B_psi(0)=0.
    dr_v <- switch(a$family,
      gaussian = 1 / unname(a$psi["sigma2"]),
      poisson = exp(eta),
      nbinom = 1 / (exp(-eta) + 1 / unname(a$psi["theta"])),
      cpg = exp((2 - unname(a$psi["power"])) * eta) /
            unname(a$psi["phi"]))
  } else if (a$family == "gaussian") {
    s2 <- unname(a$psi["sigma2"])
    dr_v <- s2 / (v + s2)^2
    dr_psi["sigma2"] <- -v / (v + s2)^2
  } else {
    # log B = 2 eta + v + log(expm1(v)).
    # Derivatives of log W are computed using stable mixture weights.
    dlogB_v <- 1 + 1 / (-expm1(-v))
    if (a$family == "poisson") {
      dlogW_eta <- 1
      dlogW_v <- 0.5
    } else if (a$family == "nbinom") {
      theta <- unname(a$psi["theta"])
      logA <- eta + v / 2
      logC <- 2 * eta + 2 * v - log(theta)
      weightC <- stats::plogis(logC - logA)
      dlogW_eta <- 1 + weightC
      dlogW_v <- 0.5 + 1.5 * weightC
      dr_psi["theta"] <- rho * (1 - rho) * weightC / theta
    } else {
      phi <- unname(a$psi["phi"])
      power <- unname(a$psi["power"])
      dlogW_eta <- power
      dlogW_v <- power^2 / 2
      dr_psi["phi"] <- -rho * (1 - rho) / phi
      dr_psi["power"] <- -rho * (1 - rho) * (eta + power * v)
    }
    dr_eta <- rho * (1 - rho) * (2 - dlogW_eta)
    dr_v <- rho * (1 - rho) * (dlogB_v - dlogW_v)
  }
  out <- c(beta0 = dr_eta, beta1 = x * dr_eta,
           dr_psi, gamma11 = dr_v, gamma12 = 2 * x * dr_v,
           gamma22 = x^2 * dr_v)
  if (any(!is.finite(out)))
    stop("Nonfinite VPC gradient; check parameter magnitudes.", call. = FALSE)
  out
}

#' Ordinary interior delta-method confidence interval for VPC or contrast
#'
#' @param fit A glmmVPC_fit object with positive-definite Gamma.
#' @param x Covariate value of rho(x), or first value for a contrast.
#' @param information Fisher information matrix in the exact parameter order
#'   beta0, beta1, psi, gamma11, gamma12, gamma22.
#' @param x1 Optional second covariate value. When supplied, target is
#'   rho(x1) - rho(x).
#' @param level Confidence level in (0,1).
#' @param information_scale "total" if information is for all clusters,
#'   "per_cluster" if it is per cluster (requires n_clusters).
#' @param n_clusters Number of independent clusters for per-cluster information.
#' @param clip Whether to intersect the interval with the target's feasible
#'   range: [0,1] for VPC, [-1,1] for contrast.
#' @return A list containing estimate, se, interval, gradient and diagnostics.
#' @export
vpc_ci_delta <- function(fit, x, information, x1 = NULL, level = 0.95,
                         information_scale = c("total", "per_cluster"),
                         n_clusters = NULL, clip = TRUE) {
  a <- .vpc_ci_params(fit)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop("`x` must be a finite scalar.", call. = FALSE)
  if (!is.null(x1) && (!is.numeric(x1) || length(x1) != 1L ||
                      !is.finite(x1) || x1 == x))
    stop("`x1` must be finite and different from x.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("`level` must lie in (0,1).", call. = FALSE)
  if (!is.logical(clip) || length(clip) != 1L || is.na(clip))
    stop("`clip` must be TRUE or FALSE.", call. = FALSE)
  eigG <- eigen(a$Gamma, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigG) <= 0)
    stop("Ordinary delta CI requires positive-definite Gamma; use constrained CI at the boundary.",
         call. = FALSE)
  information_scale <- match.arg(information_scale)
  d <- length(a$names)
  if (!is.matrix(information) || !is.numeric(information) ||
      !identical(dim(information), c(d, d)) || any(!is.finite(information)))
    stop("`information` must be a finite square matrix of the required dimension.",
         call. = FALSE)
  if (!is.null(rownames(information)) || !is.null(colnames(information))) {
    if (!identical(rownames(information), a$names) ||
        !identical(colnames(information), a$names))
      stop("Information matrix dimnames must exactly match parameter order.",
           call. = FALSE)
  } else {
    warning("Unlabeled information matrix: assuming beta, psi, gamma parameter order.",
            call. = FALSE)
  }
  if (max(abs(information - t(information))) >
      1e-8 * max(1, max(abs(information))))
    stop("Information matrix must be symmetric.", call. = FALSE)
  information <- (information + t(information)) / 2
  if (information_scale == "per_cluster") {
    if (!is.numeric(n_clusters) || length(n_clusters) != 1L ||
        !is.finite(n_clusters) || n_clusters < 1 ||
        n_clusters != floor(n_clusters))
      stop("Provide a positive integer `n_clusters`.", call. = FALSE)
    information <- information * n_clusters
  }
  cholI <- tryCatch(chol(information), error = function(e) NULL)
  if (is.null(cholI)) stop("Information matrix must be positive definite.", call. = FALSE)
  g <- vpc_gradient(fit, x)
  estimate <- as.numeric(glmm_vpc(a$beta, a$Gamma, x, a$family, a$psi))
  target <- "vpc"
  bounds <- c(0, 1)
  if (!is.null(x1)) {
    g <- vpc_gradient(fit, x1) - g
    estimate <- as.numeric(glmm_vpc(a$beta, a$Gamma, x1, a$family, a$psi)) -
      estimate
    target <- "difference"
    bounds <- c(-1, 1)
  }
  # Solve I y = g without forming I^-1.
  y <- backsolve(cholI, forwardsolve(t(cholI), unname(g)))
  variance <- sum(unname(g) * y)
  if (!is.finite(variance) || variance < -1e-10)
    stop("Invalid delta variance.", call. = FALSE)
  se <- sqrt(max(0, variance))
  z <- stats::qnorm((1 + level) / 2)
  ci <- estimate + c(-1, 1) * z * se
  if (clip) ci <- pmin(bounds[2], pmax(bounds[1], ci))
  names(ci) <- c("lower", "upper")
  structure(list(target = target, x = x, x1 = x1, estimate = estimate,
                 se = se, conf.int = ci, level = level, gradient = g,
                 variance = variance, information_scale = information_scale,
                 n_clusters = if (information_scale == "per_cluster") n_clusters else NULL,
                 method = "ordinary interior delta method",
                 boundary_valid = FALSE),
            class = "glmmVPC_ci")
}
