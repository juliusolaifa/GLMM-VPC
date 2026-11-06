# Boundary-aware, rank-adaptive constrained delta method.
# Requires vpc_ci_delta.R (vpc_gradient, .vpc_ci_params) and
# cone_geometry.R (project_psd_metric).
#
# Uses TOTAL Fisher information. Simulated displacements already have the
# estimator's scale; no additional division by sqrt(n) is performed.
# Pointwise fixed-parameter asymptotics only; no uniform coverage guarantee.

.cic_information <- function(information, names_required, scale, n_clusters) {
  d <- length(names_required)
  if (!is.matrix(information) || !is.numeric(information) ||
      !identical(dim(information), c(d, d)) || any(!is.finite(information)))
    stop("Information must be a finite square matrix in full parameter order.",
         call. = FALSE)
  if (!is.null(dimnames(information)) &&
      (!identical(rownames(information), names_required) ||
       !identical(colnames(information), names_required)))
    stop("Information dimnames must match beta, psi, gamma order.", call. = FALSE)
  if (max(abs(information - t(information))) >
      1e-8 * max(1, max(abs(information))))
    stop("Information must be symmetric.", call. = FALSE)
  I <- (information + t(information))/2
  if (scale == "per_cluster") {
    if (is.null(n_clusters)) stop("n_clusters required for per-cluster information.")
    I <- I * n_clusters
  }
  if (inherits(try(chol(I), silent = TRUE), "try-error"))
    stop("Information must be positive definite.", call. = FALSE)
  I
}

#' Select an estimated covariance-rank stratum
#'
#' For consistency at a fixed true covariance, threshold tau_n must satisfy
#' tau_n -> 0 and sqrt(n)*tau_n -> infinity. The default n^(-1/4)
#' is illustrative, not a finite-sample optimal choice.
#' @param Gamma Symmetric 2x2 covariance estimate.
#' @param n_clusters Number of independent clusters.
#' @param rank_threshold Optional nonnegative eigenvalue threshold.
#' @return List with rank, eigenvalues, null eigenvector and threshold.
#' @export
vpc_ci_covariance_rank <- function(Gamma, n_clusters, rank_threshold = NULL) {
  if (!is.matrix(Gamma) || !is.numeric(Gamma) ||
      !identical(dim(Gamma), c(2L, 2L)) || any(!is.finite(Gamma)) ||
      max(abs(Gamma-t(Gamma))) > 1e-8)
    stop("Gamma must be a finite symmetric 2x2 matrix.", call. = FALSE)
  if (!is.numeric(n_clusters) || length(n_clusters) != 1L ||
      !is.finite(n_clusters) || n_clusters < 2 ||
      n_clusters != floor(n_clusters))
    stop("n_clusters must be an integer >= 2.", call. = FALSE)
  if (is.null(rank_threshold)) rank_threshold <- n_clusters^(-1/4)
  if (!is.numeric(rank_threshold) || length(rank_threshold) != 1L ||
      !is.finite(rank_threshold) || rank_threshold < 0)
    stop("rank_threshold must be a finite nonnegative scalar.", call. = FALSE)
  ev <- eigen((Gamma+t(Gamma))/2, symmetric = TRUE)
  if (min(ev$values) < -1e-8)
    stop("Gamma is not positive semidefinite.", call. = FALSE)
  rank <- sum(ev$values > rank_threshold)
  list(rank = as.integer(rank), eigenvalues = ev$values,
       null_vector = if (rank == 1L) ev$vectors[, 2L] else NULL,
       threshold = rank_threshold)
}

#' Simulate the full Fisher-metric tangent-cone projection limit
#'
#' The nuisance coordinates are unconstrained, but correlated with gamma.
#' Covariance projection uses the efficient information Schur complement.
#' At rank one, the covariance tangent is the half-space
#' q' H_gamma q >= 0, with q the null eigenvector of Gamma.
#' @param information Positive-definite TOTAL information matrix.
#' @param rank Rank stratum (0,1,2).
#' @param null_vector Null eigenvector for rank 1.
#' @param B Number of Gaussian replicates.
#' @param seed Optional integer seed; preserves caller RNG state.
#' @return Matrix B x p of full-parameter local displacements.
#' @export
vpc_ci_projection_draws <- function(information, rank, null_vector = NULL,
                                    B = 5000L, seed = NULL) {
  if (!is.matrix(information) || !is.numeric(information) ||
      nrow(information) != ncol(information) || nrow(information) < 3L ||
      any(!is.finite(information)))
    stop("Invalid information matrix.", call. = FALSE)
  p <- nrow(information)
  if (!rank %in% 0:2 || length(rank) != 1L)
    stop("rank must be 0, 1, or 2.", call. = FALSE)
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 100L || B != floor(B))
    stop("B must be an integer >= 100.", call. = FALSE)
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
       seed != floor(seed)))
    stop("seed must be an integer.", call. = FALSE)
  I <- (information+t(information))/2
  if (max(abs(information-t(information))) > 1e-8*max(1,max(abs(I))))
    stop("Information must be symmetric.", call. = FALSE)
  U <- tryCatch(chol(I), error = function(e) NULL)
  if (is.null(U)) stop("Information must be positive definite.", call. = FALSE)
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  # Unconstrained local normal draws: covariance I^{-1}.
  uncon <- t(backsolve(U, matrix(stats::rnorm(p * B), nrow = p)))
  if (rank == 2L) return(uncon)
  ni <- seq_len(p - 3L)
  gi <- (p-2L):p
  A <- I[ni, ni, drop=FALSE]
  C <- I[ni, gi, drop=FALSE]
  J <- I[gi, gi, drop=FALSE] - crossprod(C, solve(A, C))
  J <- (J + t(J))/2
  if (rank == 1L) {
    if (!is.numeric(null_vector) || length(null_vector) != 2L ||
        any(!is.finite(null_vector)) || sum(null_vector^2) <= 0)
      stop("rank-one projection requires a null eigenvector.", call. = FALSE)
    q <- null_vector / sqrt(sum(null_vector^2))
    a <- c(q[1]^2, 2*q[1]*q[2], q[2]^2)
    Ja <- solve(J, a)
    denom <- as.numeric(crossprod(a, Ja))
  }
  out <- uncon
  for (b in seq_len(B)) {
    g <- uncon[b, gi]
    pg <- if (rank == 0L) {
      unname(project_psd_metric(g, J)$gamma)
    } else {
      violation <- sum(a*g)
      if (violation >= 0) g else g - violation * Ja / denom
    }
    out[b, gi] <- pg
    # Optimal nuisance displacement given the constrained gamma displacement.
    out[b, ni] <- uncon[b, ni] - drop(solve(A, C %*% (pg-g)))
  }
  out
}

#' Rank-adaptive constrained-delta confidence interval
#'
#' This simulates the projected-Gaussian limiting distribution for the full
#' parameter vector, then applies the VPC/contrast gradient. It is justified
#' pointwise at fixed covariance rank when rank selection and numerical
#' projections are consistent. Quantile atoms can invalidate equal-tail
#' coverage; in that case this function stops instead of reporting a
#' nominally calibrated interval.
#'
#' @param fit glmmVPC_fit.
#' @param x Covariate value (or baseline value for a contrast).
#' @param information Full Fisher information, in beta, psi, gamma order.
#' @param x1 Optional comparison covariate: rho(x1)-rho(x).
#' @param n_clusters Independent cluster count.
#' @param level Confidence level.
#' @param information_scale "total" or "per_cluster".
#' @param rank_threshold Optional covariance eigenvalue threshold.
#' @param B Number of Gaussian projection replicates.
#' @param seed Optional RNG seed.
#' @param clip Intersect with [0,1] or [-1,1].
#' @param atom_tol Threshold for detecting a zero atom.
#' @return A glmmVPC_ci object with simulation diagnostics.
#' @export
vpc_ci_constrained <- function(fit, x, information, x1 = NULL,
                               n_clusters, level = 0.95,
                               information_scale = c("total","per_cluster"),
                               rank_threshold = NULL, B = 5000L,
                               seed = NULL, clip = TRUE,
                               atom_tol = 1e-10) {
  a <- .vpc_ci_params(fit)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop("x must be finite scalar.", call. = FALSE)
  if (!is.null(x1) && (!is.numeric(x1) || length(x1) != 1L ||
                      !is.finite(x1) || x1 == x))
    stop("x1 must be finite and different from x.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1)
    stop("level must be in (0,1).", call. = FALSE)
  if (!is.logical(clip) || length(clip) != 1L || is.na(clip))
    stop("clip must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(atom_tol) || length(atom_tol) != 1L ||
      !is.finite(atom_tol) || atom_tol < 0)
    stop("atom_tol must be nonnegative.", call. = FALSE)
  stratum <- vpc_ci_covariance_rank(a$Gamma, n_clusters, rank_threshold)
  information_scale <- match.arg(information_scale)
  I <- .cic_information(information, a$names, information_scale, n_clusters)
  draws <- vpc_ci_projection_draws(I, stratum$rank, stratum$null_vector,
                                   B = B, seed = seed)
  grad <- vpc_gradient(fit, x)
  estimate <- as.numeric(glmm_vpc(a$beta, a$Gamma, x, a$family, a$psi))
  bounds <- c(0, 1)
  target <- "vpc"
  if (!is.null(x1)) {
    grad <- vpc_gradient(fit, x1) - grad
    estimate <- as.numeric(glmm_vpc(a$beta, a$Gamma, x1, a$family, a$psi)) -
      estimate
    bounds <- c(-1, 1)
    target <- "difference"
  }
  vals <- drop(draws %*% unname(grad))
  alpha <- 1 - level
  # Equal-tail quantiles are not automatically calibrated at jumps.
  # Detect jumps at zero that straddle a requested quantile.
  p_neg <- mean(vals < -atom_tol)
  p_nonpos <- mean(vals <= atom_tol)
  probs <- c(alpha/2, 1-alpha/2)
  if (any(probs > p_neg & probs < p_nonpos))
    stop("Equal-tail quantile falls in a zero atom; calibrated one-sided/acceptance-region inference is required.",
         call. = FALSE)
  q <- as.numeric(stats::quantile(vals, probs = probs, names = FALSE, type = 1))
  ci <- estimate - rev(q)
  if (clip) ci <- pmax(bounds[1], pmin(bounds[2], ci))
  names(ci) <- c("lower", "upper")
  structure(list(target=target, x=x, x1=x1, estimate=estimate,
                 conf.int=ci, level=level, gradient=grad,
                 rank=stratum$rank, rank_threshold=stratum$threshold,
                 eigenvalues=stratum$eigenvalues, B=B,
                 zero_atom_estimate=mean(abs(vals)<=atom_tol),
                 quantiles=setNames(q,c("lower","upper")),
                 information_scale=information_scale,
                 method="rank-adaptive constrained delta (pointwise asymptotic)",
                 uniform_coverage=FALSE),
            class="glmmVPC_ci")
}
