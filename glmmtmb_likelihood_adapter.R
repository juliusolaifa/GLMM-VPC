# glmmTMB likelihood and exact omnibus-null adapter for glmmVPC.
# Requires fit_glmm.R and constrained_fit.R.

.vpc_glmmtmb_assert_fit <- function(fit) {
  if (!inherits(fit, "glmmVPC_fit") || !identical(fit$engine, "glmmTMB"))
    stop("Requires a glmmVPC_fit produced by the glmmTMB engine.")
  if (is.null(fit$raw_fit) || is.null(fit$control))
    stop("Incomplete glmmVPC_fit: raw_fit and control are required.")
  invisible(TRUE)
}

.vpc_glmmtmb_ll_tol <- function(fit) {
  tol <- fit$control$ll_tol
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("fit$control$ll_tol must be finite and nonnegative.")
  tol
}

.vpc_glmmtmb_objective_check <- function(raw, ll, tol) {
  if (!is.numeric(ll) || length(ll) != 1L || !is.finite(ll))
    stop("Reported log likelihood must be finite and scalar.")
  if (is.null(raw$obj) || !is.function(raw$obj$fn) ||
      !is.function(raw$obj$gr) || is.null(raw$fit$par))
    stop("glmmTMB fit does not expose its native objective, gradient and parameters.")
  p <- raw$fit$par
  if (!is.numeric(p) || any(!is.finite(p)))
    stop("Invalid native parameter vector.")
  val <- tryCatch(as.numeric(raw$obj$fn(unname(p))), error = function(e) NA_real_)
  if (length(val) != 1L || !is.finite(val))
    stop("Native marginal objective is nonfinite at the selected fit.")
  gap <- abs(val + ll)
  threshold <- tol * (1 + abs(ll))
  list(objective = val, logLik = ll, absolute_gap = gap,
       tolerance = threshold, consistent = gap <= threshold)
}

# The comparison must use actual analysis rows, not only nobs.
# glmmTMB normally retains the model frame; if unavailable, do not certify.
.vpc_glmmtmb_rows <- function(raw) {
  fr <- raw$frame
  if (!is.data.frame(fr) || is.null(rownames(fr))) return(NULL)
  rownames(fr)
}

.vpc_glmmtmb_compare_rows <- function(full, restricted) {
  a <- .vpc_glmmtmb_rows(full)
  b <- .vpc_glmmtmb_rows(restricted)
  n1 <- tryCatch(as.integer(stats::nobs(full)), error = function(e) NA_integer_)
  n0 <- tryCatch(as.integer(stats::nobs(restricted)), error = function(e) NA_integer_)
  if (!is.finite(n1) || !is.finite(n0) || n1 != n0)
    stop("Full and restricted models use different numbers of observations.")
  if (is.null(a) || is.null(b))
    stop("Cannot certify identical observation rows: model frame unavailable.")
  if (length(a) != n1 || length(b) != n0 || !identical(a, b))
    stop("Full and restricted models do not use identical observations in the same order.")
  list(identical_rows = TRUE, nobs = n1)
}

#' Native glmmTMB marginal likelihood in interior coordinates
#' @export
vpc_glmmtmb_likelihood <- function(fit) {
  .vpc_glmmtmb_assert_fit(fit)
  raw <- fit$raw_fit
  check <- .vpc_glmmtmb_objective_check(raw, fit$logLik,
                                        .vpc_glmmtmb_ll_tol(fit))
  if (!check$consistent)
    stop("Native marginal objective and reported log likelihood disagree.")
  par0 <- raw$fit$par
  nll <- function(par) {
    if (!is.numeric(par) || length(par) != length(par0) ||
        any(!is.finite(par))) return(Inf)
    val <- tryCatch(raw$obj$fn(unname(par)), error = function(e) Inf)
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val))
      Inf else as.numeric(val)
  }
  gradient <- function(par) {
    if (!is.numeric(par) || length(par) != length(par0) ||
        any(!is.finite(par))) return(rep(NA_real_, length(par0)))
    ans <- tryCatch(as.numeric(raw$obj$gr(unname(par))),
                    error = function(e) rep(NA_real_, length(par0)))
    if (length(ans) != length(par0)) rep(NA_real_, length(par0)) else ans
  }
  list(nll = nll, gradient = gradient, start = par0,
       parameter_names = names(par0), engine = "glmmTMB",
       family = fit$family, formula = fit$formula,
       objective_at_fit = check$objective,
       objective_consistency = check,
       covariance_parameterization = "glmmTMB-native; positive definite only")
}

#' Extract family-specific parameters from a restricted glmmTMB fit
.glmmvpc_extract_null_family_params <- function(raw_fit, family) {
  pars <- switch(family,
    gaussian = c(phi = unname(stats::sigma(raw_fit)^2)),
    poisson = numeric(),
    nbinom = c(theta = unname(stats::sigma(raw_fit))),
    cpg = {
      psi <- raw_fit$fit$par[names(raw_fit$fit$par) == "psi"]
      if (length(psi) != 1L || !is.finite(psi))
        stop("Unable to extract Tweedie power from null fit.")
      c(phi = unname(stats::sigma(raw_fit)),
        power = 1 + stats::plogis(unname(psi)))
    },
    stop("Unsupported family: ", family)
  )
  if (any(!is.finite(pars)) || any(pars <= 0) ||
      (family == "cpg" && !(pars[["power"]] > 1 && pars[["power"]] < 2)))
    stop("Invalid restricted family parameters.")
  pars
}

#' Exact omnibus null fit: Gamma = 0
#' @export
fit_vpc_omnibus_null <- function(fit, require_acceptable = TRUE) {
  .vpc_glmmtmb_assert_fit(fit)
  if (!is.logical(require_acceptable) || length(require_acceptable) != 1L ||
      is.na(require_acceptable))
    stop("require_acceptable must be TRUE or FALSE.")
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Install lme4 for formula reduction.")
  if (!exists(".glmmvpc_glmmtmb_adapter", mode = "function"))
    stop("Source fit_glmm.R first.")
  null_formula <- lme4::nobars(fit$formula)
  environment(null_formula) <- environment(fit$formula)
  ans <- .glmmvpc_glmmtmb_adapter(null_formula, fit$data,
                                   fit$family, fit$control)
  rows <- .vpc_glmmtmb_compare_rows(fit$raw_fit, ans$fit)
  d <- ans$diagnostics
  ll <- as.numeric(d$logLik)
  check <- .vpc_glmmtmb_objective_check(ans$fit, ll,
                                        .vpc_glmmtmb_ll_tol(fit))
  if (!check$consistent)
    stop("Omnibus null objective and log likelihood are inconsistent.")
  if (require_acceptable && !isTRUE(d$numerically_acceptable))
    stop("Omnibus null is not numerically acceptable; inspect diagnostics.")
  beta <- glmmTMB::fixef(ans$fit)$cond
  pars <- .glmmvpc_extract_null_family_params(ans$fit, fit$family)
  d$feasible <- TRUE
  d$constraint_violation <- 0
  d$constraint <- "exact omission of random effects"
  d$observation_rows_identical <- rows$identical_rows
  d$nobs_compared <- rows$nobs
  d$native_objective_consistency <- check
  d$covariance_strata_checked <- "exact vertex: Gamma = 0"
  d$boundary_complete <- TRUE  # For this restricted null ONLY.
  d$global_optimality_verified <- FALSE
  structure(list(
    label = "omnibus: Gamma = 0", par = c(beta, pars),
    free_par = ans$fit$fit$par, beta = beta,
    Gamma = matrix(0, 2L, 2L,
                   dimnames = list(fit$random_design_columns,
                                   fit$random_design_columns)),
    family_params = pars, logLik = ll, objective = check$objective,
    formula = null_formula, engine = "glmmTMB", raw_fit = ans$fit,
    diagnostics = d, attempts = ans$attempts
  ), class = "glmmVPC_constrained_fit")
}

#' Observed omnibus LRT, without asymptotic calibration
#' @export
vpc_omnibus_lrt <- function(fit, require_acceptable = TRUE,
                            negative_tol = 1e-6) {
  .vpc_glmmtmb_assert_fit(fit)
  if (require_acceptable &&
      !isTRUE(fit$diagnostics$numerically_acceptable))
    stop("Full fit is not numerically acceptable.")
  null_fit <- fit_vpc_omnibus_null(fit, require_acceptable)
  stat <- vpc_lrt_statistic(fit, null_fit, negative_tol)
  list(statistic = stat, full_fit = fit, null_fit = null_fit,
       full_logLik = fit$logLik, null_logLik = null_fit$logLik,
       calibration = "not computed",
       diagnostics = list(full = fit$diagnostics,
                          null = null_fit$diagnostics))
}
