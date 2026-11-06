# Public confidence-interval dispatcher for glmmVPC.
# This file delegates statistical calculations to the three canonical
# method modules. No new calibration or fitting method is introduced.

#' Confidence interval or confidence set for a VPC or VPC contrast
#'
#' @param fit Fitted glmmVPC_fit object.
#' @param x Finite covariate value.
#' @param x1 Optional second covariate; target is rho(x1)-rho(x).
#' @param method One of "delta", "constrained", "profile".
#' @param level Confidence level in (0,1).
#' @param information Full Fisher information in beta, psi, gamma order.
#'   Required for delta and constrained methods.
#' @param information_scale "total" or "per_cluster".
#' @param n_clusters Independent cluster count, required for constrained
#'   intervals and per-cluster information.
#' @param clip Intersect delta-based intervals with feasible target range.
#' @param rank_threshold Optional eigenvalue threshold for constrained method.
#' @param B Number of projected Gaussian draws for constrained method.
#' @param seed Optional seed for constrained method.
#' @param atom_tol Tolerance for zero-atom detection.
#' @param grid Candidate grid for profile inversion.
#' @param constrained_fit Candidate-specific equality-constrained fitter.
#' @param critical_value Candidate-specific critical-value callback.
#' @param require_certified Require certified profile restricted fits.
#' @param negative_tol Tolerance for negative likelihood-ratio roundoff.
#' @return For delta/constrained, a glmmVPC_ci object; for profile, a
#'   grid-based confidence-set object. Profile results are NOT continuous
#'   intervals, and their coverage depends on supplied calibration.
#' @export
vpc_ci <- function(fit, x, x1 = NULL,
                   method = c("delta", "constrained", "profile"),
                   level = 0.95, information = NULL,
                   information_scale = c("total", "per_cluster"),
                   n_clusters = NULL, clip = TRUE,
                   rank_threshold = NULL, B = 5000L,
                   seed = NULL, atom_tol = 1e-10,
                   grid = NULL, constrained_fit = NULL,
                   critical_value = NULL,
                   require_certified = TRUE,
                   negative_tol = 1e-6) {
  method <- match.arg(method)
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be a finite number in (0,1).", call. = FALSE)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop("x must be a finite scalar.", call. = FALSE)
  if (!is.null(x1) && (!is.numeric(x1) || length(x1) != 1L ||
                      !is.finite(x1) || x1 == x))
    stop("x1 must be finite and different from x.", call. = FALSE)

  if (method == "profile") {
    if (is.null(grid) || !is.function(constrained_fit) ||
        !is.function(critical_value))
      stop("Profile inversion requires grid, constrained_fit, and critical_value.",
           call. = FALSE)
    return(vpc_ci_profile(
      fit = fit, x = x, x1 = x1, grid = grid,
      constrained_fit = constrained_fit,
      critical_value = critical_value, level = level,
      require_certified = require_certified,
      negative_tol = negative_tol
    ))
  }

  if (is.null(information))
    stop("A full Fisher information matrix is required for this method.",
         call. = FALSE)
  information_scale <- match.arg(information_scale)
  if (method == "delta") {
    return(vpc_ci_delta(
      fit = fit, x = x, x1 = x1, information = information,
      level = level, information_scale = information_scale,
      n_clusters = n_clusters, clip = clip
    ))
  }

  if (is.null(n_clusters))
    stop("Constrained delta inference requires n_clusters.", call. = FALSE)
  vpc_ci_constrained(
    fit = fit, x = x, x1 = x1, information = information,
    n_clusters = n_clusters, level = level,
    information_scale = information_scale,
    rank_threshold = rank_threshold, B = B, seed = seed,
    clip = clip, atom_tol = atom_tol
  )
}
