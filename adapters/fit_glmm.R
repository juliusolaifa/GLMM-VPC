# Modular single-model GLMM fitting interface for glmmVPC.
# Initial adapter: glmmTMB. Additional engines register their own adapters.

.glmmvpc_fit_engines <- new.env(parent = emptyenv())

register_glmm_engine <- function(name, adapter, overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name),
            is.function(adapter), is.logical(overwrite), length(overwrite) == 1L,
            !is.na(overwrite))
  if (exists(name, envir = .glmmvpc_fit_engines, inherits = FALSE) && !overwrite)
    stop("Engine already registered: ", name)
  assign(name, adapter, envir = .glmmvpc_fit_engines)
  invisible(name)
}

.glmmvpc_family <- function(family) {
  if (inherits(family, "family")) {
    stop("Pass a supported family name, not a family object.")
  }
  match.arg(family, c("gaussian", "poisson", "nbinom", "cpg"))
}

.glmmvpc_extract <- function(raw_fit, family, random_terms) {
  beta <- glmmTMB::fixef(raw_fit)$cond
  vc <- glmmTMB::VarCorr(raw_fit)$cond
  if (length(vc) != 1L)
    stop("Current VPC adapter requires exactly one random-effects grouping factor.")
  Gamma <- as.matrix(vc[[1L]])
  if (!identical(dim(Gamma), c(2L, 2L)))
    stop("Current VPC inference requires a 2 x 2 random-intercept-and-slope covariance matrix.")
  Gamma <- (Gamma + t(Gamma)) / 2
  if (any(!is.finite(Gamma))) stop("Nonfinite covariance estimate.")
  fp <- switch(family,
    gaussian = c(sigma2 = unname(stats::sigma(raw_fit)^2)),
    poisson = numeric(),
    nbinom = c(theta = unname(stats::sigma(raw_fit))),
    cpg = {
      psi <- raw_fit$fit$par[names(raw_fit$fit$par) == "psi"]
      if (length(psi) != 1L || !is.finite(psi))
        stop("Unable to extract Tweedie power parameter.")
      c(phi = unname(stats::sigma(raw_fit)), power = 1 + stats::plogis(unname(psi)))
    }
  )
  if (any(!is.finite(fp)) || any(fp <= 0)) stop("Invalid family parameters.")
  if (family == "cpg" && !(fp[["power"]] > 1 && fp[["power"]] < 2))
    stop("Tweedie power must lie strictly between 1 and 2.")
  list(beta = beta, Gamma = Gamma, family_params = fp)
}

.glmmvpc_diagnose <- function(fit, warnings, grad_tol, ll_tol, elapsed) {
  ll <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
  obj <- tryCatch(as.numeric(fit$fit$objective), error = function(e) NA_real_)
  gr <- tryCatch(fit$obj$gr(fit$fit$par), error = function(e) NA_real_)
  maxgrad <- if (length(gr) && all(is.finite(gr))) max(abs(gr)) else NA_real_
  optimizer_ok <- isTRUE(fit$fit$convergence == 0L)
  hessian_ok <- isTRUE(fit$sdr$pdHess)
  likelihood_ok <- length(ll) == 1L && is.finite(ll)
  objective_consistent <- length(obj) == 1L && is.finite(obj) &&
    likelihood_ok && abs(ll + obj) <= ll_tol * (1 + abs(ll))
  stationary <- is.finite(maxgrad) && maxgrad <= grad_tol
  list(convergence_code = fit$fit$convergence,
       optimizer_message = fit$fit$message,
       optimizer_ok = optimizer_ok, hessian_ok = hessian_ok,
       likelihood_ok = likelihood_ok,
       objective_consistent = objective_consistent,
       stationary = stationary, max_abs_gradient = maxgrad,
       logLik = ll, objective = obj,
       numerically_acceptable = optimizer_ok && hessian_ok &&
         objective_consistent && stationary,
       elapsed_seconds = elapsed, warnings = warnings)
}

.glmmvpc_glmmtmb_adapter <- function(formula, data, family, control) {
  if (!requireNamespace("glmmTMB", quietly = TRUE))
    stop("Install 'glmmTMB' to use this engine.")
  fam <- switch(family,
    gaussian = stats::gaussian(link = "identity"),
    poisson = stats::poisson(link = "log"),
    nbinom = glmmTMB::nbinom2(link = "log"),
    cpg = glmmTMB::tweedie(link = "log"))

  fit_one <- function(optimizer) {
    warnings <- character()
    t0 <- proc.time()[["elapsed"]]
    out <- tryCatch(withCallingHandlers({
      ctrl <- if (optimizer == "nlminb") {
        glmmTMB::glmmTMBControl(
          optimizer = stats::nlminb,
          optCtrl = list(iter.max = control$nlminb_iter_max,
                         eval.max = control$nlminb_eval_max))
      } else {
        glmmTMB::glmmTMBControl(
          optimizer = stats::optim,
          optArgs = list(method = "BFGS"),
          optCtrl = list(maxit = control$bfgs_maxit))
      }
      glmmTMB::glmmTMB(formula = formula, data = data, family = fam,
                       control = ctrl)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (inherits(out, "error"))
      return(list(fit = NULL, diagnostics = NULL,
                  error = conditionMessage(out), warnings = warnings,
                  elapsed_seconds = elapsed))
    list(fit = out,
         diagnostics = .glmmvpc_diagnose(out, warnings, control$grad_tol,
                                         control$ll_tol, elapsed), error = NULL)
  }

  candidates <- list(nlminb = fit_one("nlminb"))
  d <- candidates$nlminb$diagnostics
  retry <- is.null(d) || !isTRUE(d$optimizer_ok) ||
    !isTRUE(d$hessian_ok) || !isTRUE(d$objective_consistent) ||
    (control$retry_mode == "strict" && !isTRUE(d$stationary))
  if (control$refit && retry) candidates$bfgs <- fit_one("bfgs")

  valid <- names(candidates)[vapply(candidates, function(z)
    !is.null(z$fit) && isTRUE(z$diagnostics$likelihood_ok) &&
      isTRUE(z$diagnostics$objective_consistent), logical(1))]
  if (!length(valid)) {
    errors <- vapply(candidates, function(z) {
      if (!is.null(z$error)) z$error else
        paste0("invalid likelihood/consistency; convergence=",
               z$diagnostics$convergence_code)
    }, character(1))
    stop("No finite, objective-consistent fit: ",
         paste(names(errors), errors, sep = ": ", collapse = "; "))
  }
  # Prefer a numerically acceptable candidate. Within that class, select
  # the highest likelihood. Retain the remaining attempts for diagnostics.
  acceptable <- valid[vapply(candidates[valid], function(z)
    isTRUE(z$diagnostics$numerically_acceptable), logical(1))]
  pool <- if (length(acceptable)) acceptable else valid
  lls <- vapply(candidates[pool], function(z) z$diagnostics$logLik, numeric(1))
  selected <- pool[which.max(lls)]
  ans <- candidates[[selected]]
  ans$diagnostics$selected_optimizer <- selected
  ans$diagnostics$fit_attempts <- length(candidates)
  ans$diagnostics$unresolved_convergence <-
    !isTRUE(ans$diagnostics$numerically_acceptable)
  ans$attempts <- lapply(candidates, function(z)
    if (is.null(z$diagnostics)) list(error = z$error, warnings = z$warnings)
    else z$diagnostics)
  ans
}

register_glmm_engine("glmmTMB", .glmmvpc_glmmtmb_adapter)

#' Fit one random-intercept-and-slope GLMM
#'
#' @param formula Full mixed-model formula, e.g. y ~ x + (1 + x | id).
#' @param data Data frame.
#' @param family gaussian, poisson, nbinom (NB2), or cpg (Tweedie).
#' @param engine Registered fitting engine (currently glmmTMB).
#' @param grad_tol Maximum absolute gradient for stationarity diagnostics.
#' @param singular_tol Relative covariance eigenvalue tolerance.
#' @param refit Retry with BFGS when the initial fit has numerical problems.
#' @param retry_mode "failed" or "strict" (also retry nonstationary fits).
#' @param bfgs_maxit Maximum BFGS iterations.
#' @param nlminb_iter_max Maximum nlminb iterations.
#' @param nlminb_eval_max Maximum nlminb evaluations.
#' @param ll_tol Likelihood/objective consistency tolerance.
#' @return An object of class glmmVPC_fit.
#' @export
fit_glmm <- function(formula, data, family = c("gaussian", "poisson", "nbinom", "cpg"),
                     engine = "glmmTMB", grad_tol = 1e-4, singular_tol = 1e-6,
                     refit = TRUE, retry_mode = c("failed", "strict"),
                     bfgs_maxit = 300L, nlminb_iter_max = 2000L,
                     nlminb_eval_max = 4000L, ll_tol = 1e-6) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("formula must be a two-sided mixed-model formula.")
  if (!is.data.frame(data)) stop("data must be a data frame.")
  family <- .glmmvpc_family(family)
  retry_mode <- match.arg(retry_mode)
  if (!is.character(engine) || length(engine) != 1L ||
      !exists(engine, envir = .glmmvpc_fit_engines, inherits = FALSE))
    stop("Unknown engine: ", paste(engine, collapse = ", "))
  if (!is.logical(refit) || length(refit) != 1L || is.na(refit))
    stop("refit must be TRUE or FALSE.")
  for (nm in c("grad_tol", "singular_tol", "ll_tol")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
        v < 0 || (nm == "grad_tol" && v == 0))
      stop("Invalid ", nm)
  }
  for (nm in c("bfgs_maxit", "nlminb_iter_max", "nlminb_eval_max")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
        v < 1 || v != as.integer(v)) stop("Invalid ", nm)
  }
  # Scope restriction: exactly one 2D random-effects block.
  # Validate formula structure through lme4-style formula parser in glmmTMB.
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required to validate random-effects structure.")
  bars <- lme4::findbars(formula)
  if (length(bars) != 1L)
    stop("Exactly one random-effects term is required: (1 + x | id).")
  random_formula <- stats::as.formula(call("~", bars[[1L]][[2L]]),
                                     env = environment(formula))
  re_design <- stats::model.matrix(random_formula, data = data)
  if (ncol(re_design) != 2L || !"(Intercept)" %in% colnames(re_design))
    stop("Random effects must consist of an intercept and one slope.")
  if (anyNA(re_design) || qr(re_design)$rank < 2L)
    stop("Random-effects design must have no missing values and rank two.")
  # Reject uncorrelated || specifications: covariance must be unrestricted PSD.
  if (grepl("\\|\\|", paste(deparse(formula), collapse = "")))
    stop("Use a correlated random-effects term (1 + x | id), not ||.")
  group_name <- paste(deparse(bars[[1L]][[3L]]), collapse = "")
  if (!group_name %in% names(data))
    stop("Grouping factor must be a single column in data.")
  if (anyNA(data[[group_name]]) || length(unique(data[[group_name]])) < 2L)
    stop("Grouping factor must contain at least two observed groups and no missing values.")

  control <- list(grad_tol = grad_tol, ll_tol = ll_tol, refit = refit,
                  retry_mode = retry_mode, bfgs_maxit = as.integer(bfgs_maxit),
                  nlminb_iter_max = as.integer(nlminb_iter_max),
                  nlminb_eval_max = as.integer(nlminb_eval_max))
  adapter <- get(engine, envir = .glmmvpc_fit_engines, inherits = FALSE)
  result <- adapter(formula, data, family, control)
  if (!is.list(result) || is.null(result$fit) || is.null(result$diagnostics))
    stop("Engine adapter must return a fit and diagnostics.")
  # Custom adapters supply normalized parameters; glmmTMB uses built-in extraction.
  pars <- if (engine == "glmmTMB")
    .glmmvpc_extract(result$fit, family, bars) else result$parameters
  if (!is.list(pars) || is.null(pars$beta) || is.null(pars$Gamma) ||
      is.null(pars$family_params) || !identical(dim(pars$Gamma), c(2L, 2L)))
    stop("Adapter must supply beta, 2x2 Gamma, and family_params.")
  eig <- eigen(pars$Gamma, symmetric = TRUE, only.values = TRUE)$values
  result$diagnostics$min_eigenvalue <- min(eig)
  result$diagnostics$max_eigenvalue <- max(eig)
  result$diagnostics$near_singular <-
    min(eig) <= singular_tol * max(1, max(eig))
  # A numerically acceptable optimizer result is not a certificate that
  # all PSD covariance strata have been explored.
  result$diagnostics$boundary_complete <- FALSE
  result$diagnostics$covariance_strata_checked <- "native glmmTMB interior"
  structure(list(
    beta = pars$beta,
    Gamma = pars$Gamma,
    family = family,
    family_params = pars$family_params,
    logLik = result$diagnostics$logLik,
    diagnostics = result$diagnostics,
    refit_diagnostics = result$attempts,
    engine = engine,
    formula = formula,
    data = data,
    group = group_name,
    random_design_columns = colnames(re_design),
    raw_fit = result$fit,
    control = control
  ), class = "glmmVPC_fit")
}

#' @export
print.glmmVPC_fit <- function(x, ...) {
  cat("glmmVPC fit\n", "Family: ", x$family,
      " | Engine: ", x$engine, "\n", sep = "")
  cat("logLik: ", format(x$logLik, digits = 8),
      " | Numerically acceptable: ",
      isTRUE(x$diagnostics$numerically_acceptable), "\n", sep = "")
  invisible(x)
}
