# Difference-contrast LRT adapter for the glmmTMB engine.
# IMPORTANT: glmmTMB native theta coordinates cover Gamma > 0 only.
# This is an INTERIOR constrained optimizer, not a complete PSD-stratified
# global null fitter. It must not be used for final boundary-aware inference.
# Requires fit_glmm.R, constrained_fit.R, glmmtmb_likelihood_adapter.R.

.vpc_diff_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop(name, ' must be one finite number.')
  as.numeric(x)
}

#' Response-scale VPC from marginal mean and variance components
#' @export
vpc_response_value <- function(beta, Gamma, family, family_params, x,
                               fixed_design = NULL) {
  x <- .vpc_diff_scalar(x, 'x')
  if (!is.matrix(Gamma) || !identical(dim(Gamma), c(2L, 2L)) ||
      any(!is.finite(Gamma))) stop('Gamma must be a finite 2 x 2 matrix.')
  if (min(eigen((Gamma + t(Gamma))/2, symmetric=TRUE,
                only.values=TRUE)$values) < -1e-8)
    stop('Gamma is not positive semidefinite.')
  if (is.null(fixed_design)) {
    if (length(beta) != 2L)
      stop('For >2 fixed coefficients, provide fixed_design(x, beta).')
    eta <- unname(beta[1L] + beta[2L] * x)
  } else {
    eta <- fixed_design(x, beta)
    eta <- .vpc_diff_scalar(eta, 'fixed linear predictor')
  }
  v <- max(0, drop(crossprod(c(1, x), Gamma %*% c(1, x))))
  if (family == 'gaussian') {
    phi <- unname(family_params['phi'])
    if (length(phi) != 1L || !is.finite(phi) || phi <= 0)
      stop('Gaussian residual variance phi must be positive.')
    return(v/(v+phi))
  }
  # Log-link families. Scale by exp(2*eta+v) to reduce overflow.
  between <- -expm1(-v)  # B / exp(2 eta + 2 v)
  within <- switch(family,
                   poisson = exp(-eta - 1.5*v),
                   nbinom = {
                     theta <- unname(family_params['theta'])
                     if (length(theta) != 1L || !is.finite(theta) || theta <= 0)
                       stop('NB2 theta must be positive.')
                     exp(-eta - 1.5*v) + 1/theta
                   },
                   cpg = {
                     phi <- unname(family_params['phi'])
                     power <- unname(family_params['power'])
                     if (length(phi) != 1L || !is.finite(phi) || phi <= 0 ||
                         length(power) != 1L || !is.finite(power) ||
                         power <= 1 || power >= 2) stop('Invalid CPG phi/power.')
                     phi*exp((power-2)*eta + (0.5*power^2-2)*v)
                   }, stop('Unsupported family: ', family))
  between/(between+within)
}

#' Response-scale VPC difference rho(x1) - rho(x0)
#' @export
vpc_response_difference <- function(beta, Gamma, family, family_params,
                                    x0, x1, fixed_design = NULL) {
  vpc_response_value(beta, Gamma, family, family_params, x1, fixed_design) -
    vpc_response_value(beta, Gamma, family, family_params, x0, fixed_design)
}

# Extract native glmmTMB coordinates for a single unstructured 2D RE term.
.vpc_diff_unpack <- function(par, fit) {
  nm <- names(fit$raw_fit$fit$par)
  if (length(par) != length(nm)) stop('Parameter length mismatch.')
  beta <- par[nm == 'beta']
  names(beta) <- names(fit$beta)
  if (length(beta) != length(fit$beta)) stop('Unexpected beta layout.')
  theta <- par[nm == 'theta']
  if (length(theta) != 3L)
    stop('Expected three glmmTMB unstructured random-effect parameters.')
  sd1 <- exp(theta[1]); sd2 <- exp(theta[2])
  corr <- theta[3]/sqrt(1+theta[3]^2)
  Gamma <- matrix(c(sd1^2, corr*sd1*sd2,
                    corr*sd1*sd2, sd2^2), 2, 2)
  disp <- par[nm == 'betadisp']
  psi <- par[nm == 'psi']
  fp <- switch(fit$family,
               gaussian = {
                 if (length(disp)!=1L) stop('Expected one Gaussian dispersion.')
                 c(phi=exp(disp))
               },
               poisson = numeric(),
               nbinom = {
                 if (length(disp)!=1L) stop('Expected one NB2 dispersion.')
                 c(theta=exp(disp))
               },
               cpg = {
                 if (length(disp)!=1L || length(psi)!=1L)
                   stop('Expected Tweedie dispersion and power parameters.')
                 c(phi=exp(disp), power=1+stats::plogis(psi))
               })
  list(beta=beta, Gamma=Gamma, family_params=fp)
}

#' Fit a nonlinear VPC contrast constraint in native glmmTMB coordinates
#'
#' This is a provisional INTERIOR fit: the optimizer cannot attain rank-one
#' or vertex covariance matrices in finite theta coordinates. For the full
#' theoretical null, a PSD-stratified likelihood engine is still required.
#' @param alternative two.sided, greater or less.
#' @param fixed_design Optional function(x,beta) computing eta(x).
#' @param feasibility_tol Numerical tolerance for constraint satisfaction.
#' @export
fit_vpc_difference_interior <- function(fit, x0, x1,
                                    alternative=c('two.sided','greater','less'),
                                    fixed_design=NULL,
                                    feasibility_tol=1e-7,
                                    maxeval=3000L,
                                    require_acceptable=TRUE) {
  alternative <- match.arg(alternative)
  x0 <- .vpc_diff_scalar(x0, 'x0'); x1 <- .vpc_diff_scalar(x1, 'x1')
  if (x0 == x1) stop('x0 and x1 must be distinct.')
  if (!inherits(fit,'glmmVPC_fit') || !identical(fit$engine,'glmmTMB'))
    stop('Requires a glmmTMB glmmVPC_fit.')
  if (!requireNamespace('nloptr', quietly=TRUE))
    stop('Install nloptr for constrained difference fitting.')
  if (!is.null(fixed_design) && !is.function(fixed_design))
    stop('fixed_design must be a function(x,beta).')
  if (require_acceptable && !isTRUE(fit$diagnostics$numerically_acceptable))
    stop('Full fit has unacceptable numerical diagnostics.')
  native <- vpc_glmmtmb_likelihood(fit)
  start <- unname(native$start)
  contrast <- function(par) {
    q <- .vpc_diff_unpack(par,fit)
    vpc_response_difference(q$beta,q$Gamma,fit$family,q$family_params,
                            x0,x1,fixed_design)
  }
  # SLSQP uses h(x)<=0 for inequalities, h(x)=0 for equality.
  constraint <- function(par) {
    d <- contrast(par)
    switch(alternative,
           two.sided=d,
           greater=d, # H0: D <= 0
           less=-d)   # H0: D >= 0
  }
  constraint_jac <- function(par) {
    vapply(seq_along(par), function(k) {
      h <- 1e-5 * max(1,abs(par[k]))
      plus <- minus <- par
      plus[k] <- plus[k] + h
      minus[k] <- minus[k] - h
      (constraint(plus)-constraint(minus))/(2*h)
    }, numeric(1))
  }
  opts <- list(algorithm='NLOPT_LD_SLSQP',
               xtol_rel=1e-9, maxeval=as.integer(maxeval),
               print_level=0)
  # Multiple starts are useful for nonconvex nonlinear constraints.
  starts <- list(start)
  # Only include original fit as a valid restricted candidate when feasible.
  attempts <- lapply(seq_along(starts), function(i) {
    tryCatch(nloptr::nloptr(
      x0=starts[[i]], eval_f=native$nll,
      eval_grad_f=native$gradient,
      eval_g_eq=if (alternative=='two.sided') constraint else NULL,
      eval_jac_g_eq=if (alternative=='two.sided') constraint_jac else NULL,
      eval_g_ineq=if (alternative!='two.sided') constraint else NULL,
      eval_jac_g_ineq=if (alternative!='two.sided') constraint_jac else NULL,
      opts=opts), error=function(e) list(error=conditionMessage(e)))
  })
  feasible <- vapply(attempts, function(a) {
    if (is.null(a$solution) || any(!is.finite(a$solution))) return(FALSE)
    z <- tryCatch(constraint(a$solution),error=function(e) Inf)
    is.finite(z) && if (alternative=='two.sided')
      abs(z)<=feasibility_tol else z<=feasibility_tol
  }, logical(1))
  if (!any(feasible))
    stop('No feasible interior constrained solution; this does NOT rule out a PSD-boundary null optimum.')
  inds <- which(feasible)
  j <- inds[which.min(vapply(attempts[inds], function(a) native$nll(a$solution),numeric(1)))]
  a <- attempts[[j]]; pars <- .vpc_diff_unpack(a$solution,fit)
  d <- contrast(a$solution)
  status_ok <- !is.null(a$status) && a$status > 0
  diagnostics <- list(
    feasible=TRUE, constraint_value=d,
    constraint_violation=if (alternative=='two.sided') abs(d) else
      if (alternative=='greater') max(0,d) else max(0,-d),
    optimizer='nloptr SLSQP', optimizer_status=a$status,
    optimizer_message=a$message, optimizer_ok=status_ok,
    numerically_acceptable=status_ok && feasible[j],
    global_optimality_verified=FALSE,
    covariance_strata_checked='positive definite only',
    boundary_complete=FALSE,
    warning='Rank-one and zero-covariance null strata are not optimized')
  if (require_acceptable && !diagnostics$numerically_acceptable)
    stop('Restricted interior optimizer did not converge acceptably.')
  structure(list(label=paste('difference',alternative), alternative=alternative,
                 x0=x0,x1=x1, contrast=d, beta=pars$beta, Gamma=pars$Gamma,
                 family_params=pars$family_params, logLik=-native$nll(a$solution),
                 objective=native$nll(a$solution), free_par=a$solution,
                 engine='glmmTMB', diagnostics=diagnostics, attempts=attempts,
                 full_parameterization='native glmmTMB interior coordinates'),
            class='glmmVPC_constrained_fit')
}


#' Fit available strata of the response-scale VPC difference null
#'
#' This is an intermediate integration, NOT a boundary-complete fitter:
#' the rank-one stratum is not optimized. The result is deliberately
#' marked unsuitable for final calibrated inference.
#' @export
fit_vpc_difference_null <- function(
    fit, x0, x1,
    alternative = c("two.sided", "greater", "less"),
    fixed_design = NULL, feasibility_tol = 1e-7,
    maxeval = 3000L, require_acceptable = TRUE) {
  alternative <- match.arg(alternative)
  x0 <- .vpc_diff_scalar(x0, "x0")
  x1 <- .vpc_diff_scalar(x1, "x1")
  if (x0 == x1) stop("x0 and x1 must be distinct.")
  if (!inherits(fit, "glmmVPC_fit") ||
      !identical(fit$engine, "glmmTMB"))
    stop("Requires a glmmTMB glmmVPC_fit.")
  if (require_acceptable &&
      !isTRUE(fit$diagnostics$numerically_acceptable))
    stop("Full fit has unacceptable numerical diagnostics.")

  # For directional inequalities, the unrestricted MLE is a restricted
  # maximizer whenever it satisfies H0, provided the full fit is valid.
  estimated <- vpc_response_difference(
    fit$beta, fit$Gamma, fit$family, fit$family_params,
    x0, x1, fixed_design)
  full_feasible <- switch(alternative,
    two.sided = abs(estimated) <= feasibility_tol,
    greater = estimated <= feasibility_tol,
    less = estimated >= -feasibility_tol)
  if (alternative != "two.sided" && full_feasible) {
    d <- fit$diagnostics
    d$feasible <- TRUE
    d$constraint_violation <- 0
    d$selected_stratum <- "unrestricted full fit satisfying inequality"
    d$covariance_strata_checked <- d$selected_stratum
    d$boundary_complete <- TRUE
    d$global_optimality_verified <- FALSE
    d$warning <- "LRT is zero if unrestricted fit is the numerical global MLE."
    return(structure(list(
      label = paste("difference", alternative), alternative = alternative,
      x0 = x0, x1 = x1, contrast = estimated,
      beta = fit$beta, Gamma = fit$Gamma,
      family_params = fit$family_params, logLik = fit$logLik,
      objective = -fit$logLik, engine = fit$engine,
      diagnostics = d, full_fit = fit
    ), class = "glmmVPC_constrained_fit"))
  }

  # The exact vertex belongs to every difference null, including
  # both directional inequalities.
  vertex <- tryCatch(
    fit_vpc_omnibus_null(fit, require_acceptable = require_acceptable),
    error = function(e) e)
  interior <- tryCatch(
    fit_vpc_difference_interior(
      fit, x0, x1, alternative = alternative,
      fixed_design = fixed_design, feasibility_tol = feasibility_tol,
      maxeval = maxeval, require_acceptable = require_acceptable),
    error = function(e) e)

  valid <- function(z) {
    inherits(z, "glmmVPC_constrained_fit") &&
      is.numeric(z$logLik) && length(z$logLik) == 1L &&
      is.finite(z$logLik) &&
      isTRUE(z$diagnostics$feasible) &&
      (!require_acceptable ||
       isTRUE(z$diagnostics$numerically_acceptable))
  }
  ok <- c(interior = valid(interior), vertex = valid(vertex))
  if (!any(ok)) {
    msg <- function(z) if (inherits(z, "error")) conditionMessage(z) else
      "candidate failed validation"
    stop("No usable restricted candidate: interior: ", msg(interior),
         "; vertex: ", msg(vertex))
  }
  candidates <- list(interior = interior, vertex = vertex)
  lls <- vapply(names(candidates), function(k)
    if (ok[[k]]) candidates[[k]]$logLik else NA_real_, numeric(1))
  chosen <- names(which.max(replace(lls, is.na(lls), -Inf)))
  out <- candidates[[chosen]]
  out$diagnostics$selected_stratum <- chosen
  out$diagnostics$stratum_logLik <- lls
  out$diagnostics$covariance_strata_checked <- names(ok)[ok]
  out$diagnostics$stratum_errors <- lapply(candidates, function(z)
    if (inherits(z, "error")) conditionMessage(z) else NULL)
  out$diagnostics$boundary_complete <- FALSE
  out$diagnostics$global_optimality_verified <- FALSE
  out$diagnostics$warning <- paste(
    "Rank-one constrained likelihood is not optimized;",
    "this is not a certified maximum over the full PSD null.")
  out
}

#' Observed difference LRT, provisional until rank-one fitting is integrated
#' @export
vpc_difference_lrt <- function(
    fit, x0, x1,
    alternative = c("two.sided", "greater", "less"),
    fixed_design = NULL, require_acceptable = TRUE, ...) {
  alternative <- match.arg(alternative)
  null <- fit_vpc_difference_null(
    fit, x0, x1, alternative = alternative,
    fixed_design = fixed_design,
    require_acceptable = require_acceptable, ...)
  stat <- vpc_lrt_statistic(fit, null)
  list(statistic = stat, full_fit = fit, null_fit = null,
       x0 = x0, x1 = x1, alternative = alternative,
       estimated_contrast = vpc_response_difference(
         fit$beta, fit$Gamma, fit$family, fit$family_params,
         x0, x1, fixed_design),
       diagnostics = list(full = fit$diagnostics,
                          null = null$diagnostics),
       calibration = "not computed; rank-one null fitting pending")
}
