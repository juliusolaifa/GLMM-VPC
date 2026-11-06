# Profile-likelihood inversion for VPCs and VPC contrasts.
# Requires a candidate-value equality-constrained likelihood fitter and
# candidate-specific critical-value calibration supplied by the caller.
# The existing zero-null adapters DO NOT satisfy this contract for u != 0.
# No universal chi-square cutoff is assumed.

#' Profile-likelihood confidence set on a user-supplied candidate grid
#'
#' The fitter must maximize the *same marginal log-likelihood* as the full
#' fit, subject to rho(x)=u or rho(x1)-rho(x)=u, across all feasible PSD
#' covariance strata. Its returned `certified` flag must truthfully indicate
#' whether the required restricted maximum is reliable. Critical values
#' must be justified for each composite candidate null.
#'
#' @param fit A fitted glmmVPC_fit with finite `logLik`.
#' @param x Covariate value (or baseline covariate for a contrast).
#' @param x1 Optional second covariate, targeting rho(x1)-rho(x).
#' @param grid Finite sorted unique candidate values in [0,1] or [-1,1].
#' @param constrained_fit Function (fit, u, x, x1) returning a list with
#'   `logLik` (scalar), `certified` (logical scalar), and optional diagnostics.
#' @param critical_value Function (u, fit, x, x1, level) returning a finite
#'   nonnegative critical value for the candidate null.
#' @param level Confidence level.
#' @param require_certified If TRUE, abort unless each fit is certified.
#' @param negative_tol Allowed negative LRT discrepancy due to numerical error.
#' @return A grid-based confidence-set object; NOT an exact continuous
#'   interval. Gaps between grid points are not resolved.
#' @export
vpc_ci_profile <- function(fit, x, grid, constrained_fit,
                           critical_value, x1 = NULL, level = 0.95,
                           require_certified = TRUE,
                           negative_tol = 1e-6) {
  if (!inherits(fit, 'glmmVPC_fit') ||
      !is.numeric(fit$logLik) || length(fit$logLik) != 1L ||
      !is.finite(fit$logLik))
    stop('fit must be a glmmVPC_fit with finite scalar logLik.', call. = FALSE)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop('x must be a finite scalar.', call. = FALSE)
  if (!is.null(x1) && (!is.numeric(x1) || length(x1) != 1L ||
                       !is.finite(x1) || x1 == x))
    stop('x1 must be finite and different from x.', call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1)
    stop('level must be in (0,1).', call. = FALSE)
  if (!is.numeric(grid) || !length(grid) || any(!is.finite(grid)) ||
      is.unsorted(grid, strictly = TRUE))
    stop('grid must be a nonempty, strictly increasing numeric vector.',
         call. = FALSE)
  bounds <- if (is.null(x1)) c(0, 1) else c(-1, 1)
  if (any(grid < bounds[1] | grid > bounds[2]))
    stop('grid contains values outside the feasible target range.',
         call. = FALSE)
  if (!is.function(constrained_fit) || !is.function(critical_value))
    stop('Both constrained_fit and critical_value must be functions.',
         call. = FALSE)
  if (!is.logical(require_certified) || length(require_certified) != 1L ||
      is.na(require_certified))
    stop('require_certified must be TRUE or FALSE.', call. = FALSE)
  if (!is.numeric(negative_tol) || length(negative_tol) != 1L ||
      !is.finite(negative_tol) || negative_tol < 0)
    stop('negative_tol must be a finite nonnegative scalar.', call. = FALSE)

  out <- vector('list', length(grid))
  diagnostics <- vector('list', length(grid))
  for (i in seq_along(grid)) {
    u <- grid[i]
    cf <- constrained_fit(fit = fit, u = u, x = x, x1 = x1)
    if (!is.list(cf) || !is.numeric(cf$logLik) ||
        length(cf$logLik) != 1L || !is.finite(cf$logLik))
      stop(sprintf('Invalid restricted logLik at u=%.12g.', u), call. = FALSE)
    certified <- isTRUE(cf$certified)
    if (require_certified && !certified)
      stop(sprintf('Restricted maximum at u=%.12g is not certified.', u),
           call. = FALSE)
    stat <- 2 * (fit$logLik - cf$logLik)
    if (stat < -negative_tol)
      stop(sprintf('Restricted logLik exceeds full logLik at u=%.12g.', u),
           call. = FALSE)
    stat <- max(0, stat)  # Only negligible negative roundoff is truncated.
    cutoff <- critical_value(u = u, fit = fit, x = x,
                             x1 = x1, level = level)
    if (!is.numeric(cutoff) || length(cutoff) != 1L ||
        !is.finite(cutoff) || cutoff < 0)
      stop(sprintf('Invalid critical value at u=%.12g.', u), call. = FALSE)
    out[[i]] <- c(u = u, statistic = stat, critical = cutoff,
                  accepted = as.numeric(stat <= cutoff),
                  certified = as.numeric(certified))
    diagnostics[[i]] <- cf$diagnostics
  }
  tab <- as.data.frame(do.call(rbind, out))
  tab$accepted <- as.logical(tab$accepted)
  tab$certified <- as.logical(tab$certified)
  # Separate accepted grid runs; no interpolation across untested u values.
  idx <- which(tab$accepted)
  components <- list()
  if (length(idx)) {
    group <- cumsum(c(TRUE, diff(idx) != 1L))
    components <- lapply(split(idx, group), function(ii)
      c(lower = min(grid[ii]), upper = max(grid[ii])))
  }
  structure(list(target = if (is.null(x1)) 'vpc' else 'difference',
                 x = x, x1 = x1, level = level,
                 grid_results = tab, accepted_grid = grid[idx],
                 components = components, diagnostics = diagnostics,
                 all_certified = all(tab$certified),
                 method = 'grid-based profile-likelihood test inversion',
                 grid_only = TRUE,
                 coverage_guarantee = 'depends on supplied critical values and certified constrained fits'),
            class = 'glmmVPC_profile_ci')
}
