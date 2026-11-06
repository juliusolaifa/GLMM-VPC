# Pointwise zero-VPC null adapter for glmmVPC (glmmTMB engine).
# Requires fit_glmm.R, constrained_fit.R and glmmtmb_likelihood_adapter.R.
# Under Gamma >= 0, z(x0)' Gamma z(x0)=0 iff Gamma z(x0)=0.
# For random-effect design (1,x), null random effects equal u*(x-x0).
# This model covers positive-variance points on the null ray, but not its
# zero-variance vertex. An exact omnibus-null fit must also be compared.

#' Fit the pointwise zero-VPC ray away from its vertex
#' @param fit glmmVPC_fit from glmmTMB with one random intercept and slope.
#' @param x0 Finite covariate value at which VPC is tested.
#' @param require_acceptable Require numerically acceptable restricted fit.
#' @return glmmVPC_constrained_fit.
#' @export
fit_vpc_pointwise_ray <- function(fit, x0, require_acceptable = TRUE) {
  if (!inherits(fit, 'glmmVPC_fit') || !identical(fit$engine, 'glmmTMB'))
    stop('Pointwise adapter requires a glmmTMB glmmVPC_fit.')
  if (!is.numeric(x0) || length(x0) != 1L || !is.finite(x0))
    stop('x0 must be one finite numeric value.')
  if (!is.logical(require_acceptable) || length(require_acceptable) != 1L || is.na(require_acceptable))
    stop('require_acceptable must be TRUE or FALSE.')
  if (!requireNamespace('lme4', quietly = TRUE)) stop('Install lme4.')
  if (!exists('.glmmvpc_glmmtmb_adapter', mode = 'function'))
    stop('Source fit_glmm.R first.')
  bars <- lme4::findbars(fit$formula)
  if (length(bars) != 1L) stop('Exactly one random-effects term is required.')
  re_formula <- stats::as.formula(call('~', bars[[1L]][[2L]]),
                                   env = environment(fit$formula))
  Xr <- stats::model.matrix(re_formula, data = fit$data)
  if (ncol(Xr) != 2L || !identical(colnames(Xr)[1L], '(Intercept)'))
    stop('Expected random design with intercept in first column and one slope in second.')
  # This adapter uses the actual second random-design column, so it supports
  # a transformed numerical slope expression, not just a variable named x.
  slope <- Xr[, 2L]
  if (!is.numeric(slope) || any(!is.finite(slope)))
    stop('Random slope design must be finite and numeric.')
  dat <- fit$data
  temp <- '.glmmvpc_pointwise_slope'
  if (temp %in% names(dat)) stop('Reserved column name already exists: ', temp)
  dat[[temp]] <- slope - x0
  fixed <- lme4::nobars(fit$formula)
  group <- fit$group
  if (!is.character(group) || length(group) != 1L ||
      !grepl('^[A-Za-z.][A-Za-z0-9._]*$', group))
    stop('Grouping factor must be a simple syntactic column name.')
  re_term <- call('|', call('+', 0, as.name(temp)), as.name(group))
  rhs <- call('+', fixed[[3L]], call('(', re_term))
  null_formula <- stats::as.formula(call('~', fixed[[2L]], rhs),
                                     env = environment(fit$formula))
  ans <- .glmmvpc_glmmtmb_adapter(null_formula, dat, fit$family, fit$control)
  # Identical counts alone do not establish identical likelihood samples.
  rows <- .vpc_glmmtmb_compare_rows(fit$raw_fit, ans$fit)
  d <- ans$diagnostics
  if (require_acceptable && !isTRUE(d$numerically_acceptable))
    stop('Pointwise null fit is not numerically acceptable; inspect diagnostics.')
  ll <- as.numeric(d$logLik)
  if (length(ll) != 1L || !is.finite(ll))
    stop('Nonfinite restricted log likelihood.')
  check <- .vpc_glmmtmb_objective_check(ans$fit, ll,
                                        .vpc_glmmtmb_ll_tol(fit))
  if (!check$consistent)
    stop('Pointwise null objective and log likelihood disagree.')
  v <- as.matrix(glmmTMB::VarCorr(ans$fit)$cond[[group]])
  if (!identical(dim(v), c(1L, 1L)) || !is.finite(v[1L,1L]) || v[1L,1L] < 0)
    stop('Expected a finite scalar random-effect variance.')
  direction <- c(-x0, 1)
  Gamma <- unname(v[1L,1L]) * tcrossprod(direction)
  if (abs(drop(crossprod(c(1, x0), Gamma %*% c(1, x0)))) >
      1e-8 * max(1, max(abs(Gamma))))
    stop('Pointwise covariance mapping failed its independent null check.')
  dimnames(Gamma) <- list(fit$random_design_columns, fit$random_design_columns)
  beta <- glmmTMB::fixef(ans$fit)$cond
  pars <- .glmmvpc_extract_null_family_params(ans$fit, fit$family)
  d$feasible <- TRUE
  d$constraint_violation <- 0
  d$observation_rows_identical <- rows$identical_rows
  d$nobs_compared <- rows$nobs
  d$native_objective_consistency <- check
  d$covariance_strata_checked <- 'pointwise rank-one ray, positive variance only'
  # The full pointwise null is a closed ray including Gamma = 0. Native
  # glmmTMB log-SD coordinates cannot attain its endpoint exactly.
  d$boundary_complete <- FALSE
  d$global_optimality_verified <- FALSE
  d$warning <- paste('Vertex Gamma = 0 is not searched by this fit;',
                     'compare with exact omnibus null before inference.')
  structure(list(
    label = sprintf('pointwise: VPC(%s) = 0', format(x0)),
    x0 = x0, beta = beta, Gamma = Gamma, family_params = pars,
    logLik = ll, objective = check$objective, formula = null_formula,
    engine = 'glmmTMB', raw_fit = ans$fit,
    free_par = ans$fit$fit$par,
    diagnostics = within(d, {
      constraint <- 'pointwise ray: b(x) = u*(x-x0); vertex excluded'
    }),
    attempts = ans$attempts
  ), class = 'glmmVPC_constrained_fit')
}

#' Fit the closed pointwise null ray, comparing positive-variance and vertex fits
#' @description Compares the rank-one restricted fit with the exact Gamma=0
#'   fit. A failed ray optimization cannot certify that the vertex is optimal.
#' @export
fit_vpc_pointwise_null <- function(fit, x0, require_acceptable = TRUE,
                                   likelihood_tol = 1e-7) {
  .vpc_glmmtmb_assert_fit(fit)
  if (!is.numeric(x0) || length(x0) != 1L || !is.finite(x0))
    stop('x0 must be one finite numeric value.')
  if (!is.logical(require_acceptable) || length(require_acceptable) != 1L ||
      is.na(require_acceptable)) stop('require_acceptable must be TRUE or FALSE.')
  if (!is.numeric(likelihood_tol) || length(likelihood_tol) != 1L ||
      !is.finite(likelihood_tol) || likelihood_tol < 0)
    stop('likelihood_tol must be finite and nonnegative.')
  ray <- tryCatch(fit_vpc_pointwise_ray(fit, x0, require_acceptable),
                  error = function(e) e)
  vertex <- tryCatch(fit_vpc_omnibus_null(fit, require_acceptable),
                     error = function(e) e)
  good <- function(z) inherits(z, 'glmmVPC_constrained_fit') &&
    is.numeric(z$logLik) && length(z$logLik) == 1L && is.finite(z$logLik) &&
    (!require_acceptable || isTRUE(z$diagnostics$numerically_acceptable))
  ray_ok <- good(ray); vertex_ok <- good(vertex)
  if (!ray_ok && !vertex_ok)
    stop('Neither pointwise ray nor vertex fit is usable: ray: ',
         if (inherits(ray,'error')) conditionMessage(ray) else 'invalid fit',
         '; vertex: ', if (inherits(vertex,'error')) conditionMessage(vertex) else 'invalid fit')
  if (ray_ok && vertex_ok) {
    selected <- if (vertex$logLik > ray$logLik + likelihood_tol) 'vertex' else 'ray'
  } else selected <- if (ray_ok) 'ray' else 'vertex'
  ans <- if (selected == 'ray') ray else vertex
  ans$label <- sprintf('pointwise: VPC(%s) = 0', format(x0))
  ans$x0 <- x0
  ans$diagnostics$selected_stratum <- selected
  ans$diagnostics$covariance_strata_checked <- c(
    if (ray_ok) 'pointwise rank-one ray (positive variance)',
    if (vertex_ok) 'vertex (Gamma=0)')
  ans$diagnostics$boundary_complete <- ray_ok && vertex_ok
  ans$diagnostics$global_optimality_verified <- FALSE
  ans$diagnostics$stratum_logLik <- c(
    ray = if (ray_ok) ray$logLik else NA_real_,
    vertex = if (vertex_ok) vertex$logLik else NA_real_)
  ans$diagnostics$stratum_errors <- list(
    ray = if (inherits(ray,'error')) conditionMessage(ray) else NULL,
    vertex = if (inherits(vertex,'error')) conditionMessage(vertex) else NULL)
  ans$diagnostics$warning <- if (ray_ok && vertex_ok)
    'Both null strata searched; numerical convergence does not certify a global optimum.'
  else 'One null stratum was not successfully optimized; restricted maximum is not certified.'
  if (require_acceptable && !(ray_ok && vertex_ok))
    stop('Pointwise null is not boundary-complete: ', ans$diagnostics$warning)
  ans
}

