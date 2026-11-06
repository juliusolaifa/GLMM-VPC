# Generic constrained marginal-likelihood fitting for glmmVPC.
# A parameter map is exact only if its implementation is correct; independent
# validation is necessary before relying on a restricted fit for inference.

#' Define a constrained marginal-likelihood problem
#' @param nll Scalar negative log-likelihood in full parameter coordinates.
#' @param map Maps reduced coordinates to full coordinates.
#' @param start Finite reduced-coordinate starting vector.
#' @param validate Optional independent null-space validator. May return a
#'   scalar logical or a nonnegative numeric constraint violation.
#' @param label Descriptive null label.
#' @param full_names Optional full-parameter names.
#' @export
vpc_constraint <- function(nll, map, start, validate = NULL,
                           label = "constrained", full_names = NULL) {
  if (!is.function(nll) || !is.function(map)) stop("nll and map must be functions.")
  if (!is.numeric(start) || any(!is.finite(start)))
    stop("start must be a finite numeric vector.")
  if (!is.null(validate) && !is.function(validate))
    stop("validate must be NULL or a function.")
  if (!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label))
    stop("label must be a nonempty string.")
  if (!is.null(full_names) && (!is.character(full_names) || anyNA(full_names)))
    stop("full_names must be NULL or a character vector.")
  structure(list(nll=nll, map=map, start=start, validate=validate,
                 label=label, full_names=full_names),
            class="glmmVPC_constraint")
}

.vpc_constraint_violation <- function(validate, par) {
  if (is.null(validate)) return(NA_real_)
  z <- tryCatch(validate(par), error=function(e) NULL)
  if (is.logical(z) && length(z)==1L && !is.na(z))
    return(if (z) 0 else Inf)
  if (is.numeric(z) && length(z)==1L && !is.na(z) && z>=0)
    return(as.numeric(z))
  Inf
}

.vpc_reduced_fd_gradient <- function(fn, u) {
  if (!length(u)) return(numeric())
  vapply(seq_along(u), function(i) {
    h <- 1e-5*max(1,abs(u[i]))
    a <- b <- u; a[i] <- a[i]+h; b[i] <- b[i]-h
    fa <- fn(a); fb <- fn(b)
    if (!is.finite(fa) || !is.finite(fb)) return(NA_real_)
    (fa-fb)/(2*h)
  }, numeric(1))
}

#' Fit a constrained marginal likelihood
#' @param problem A vpc_constraint object.
#' @param starts Optional list of additional reduced-coordinate starts.
#' @param methods Subset of nlminb and BFGS.
#' @param grad Optional gradient of full-coordinate nll; transformed to reduced
#'   coordinates using a numerical Jacobian of map.
#' @param grad_tol Absolute reduced-gradient threshold.
#' @param maxit Iteration/evaluation budget per attempt.
#' @param feasibility_tol Maximum allowed numeric validation violation.
#' @param require_validation If TRUE, an independent validate function is required.
#' @param require_acceptable If TRUE, stop when no numerically acceptable
#'   candidate exists; otherwise return a diagnostic candidate.
#' @export
fit_glmm_constrained <- function(problem, starts=NULL,
                                 methods=c("nlminb","BFGS"),
                                 grad=NULL, grad_tol=1e-4,
                                 maxit=2000L, feasibility_tol=1e-7,
                                 require_validation=FALSE,
                                 require_acceptable=FALSE) {
  if (!inherits(problem,"glmmVPC_constraint"))
    stop("problem must be created by vpc_constraint().")
  if (!is.character(methods) || !length(methods) || anyDuplicated(methods) ||
      !all(methods %in% c("nlminb","BFGS")))
    stop("methods must be distinct selections of nlminb and/or BFGS.")
  scalar <- function(x, positive=FALSE) is.numeric(x) && length(x)==1L &&
    is.finite(x) && (if (positive) x>0 else x>=0)
  if (!scalar(grad_tol,TRUE)) stop("grad_tol must be positive and finite.")
  if (!scalar(feasibility_tol)) stop("feasibility_tol must be nonnegative and finite.")
  if (!scalar(maxit,TRUE) || maxit != floor(maxit)) stop("maxit must be a positive integer.")
  if (!is.null(grad) && !is.function(grad)) stop("grad must be a function.")
  for (v in list(require_validation,require_acceptable))
    if (!is.logical(v) || length(v)!=1L || is.na(v))
      stop("require_validation and require_acceptable must be scalar logicals.")
  if (require_validation && is.null(problem$validate))
    stop("Independent validate() is required for certified feasibility.")
  if (!is.null(starts) && !is.list(starts)) stop("starts must be a list.")
  all_starts <- c(list(problem$start),starts)
  if (!all(vapply(all_starts, function(s) is.numeric(s) &&
      length(s)==length(problem$start) && all(is.finite(s)),logical(1))))
    stop("All starts must be finite and have the same length.")

  evaluate <- function(u) {
    par <- tryCatch(problem$map(u),error=function(e) NULL)
    if (!is.numeric(par) || !length(par) || any(!is.finite(par)))
      return(list(objective=Inf,violation=Inf,par=NULL))
    if (!is.null(problem$full_names) && length(par)!=length(problem$full_names))
      return(list(objective=Inf,violation=Inf,par=NULL))
    violation <- .vpc_constraint_violation(problem$validate,par)
    if (!is.na(violation) && violation>feasibility_tol)
      return(list(objective=Inf,violation=violation,par=par))
    value <- tryCatch(problem$nll(par),error=function(e) Inf)
    if (!is.numeric(value) || length(value)!=1L || !is.finite(value)) value <- Inf
    list(objective=as.numeric(value),violation=violation,par=par)
  }
  reduced <- function(u) evaluate(u)$objective
  # Only use full-coordinate gradient for diagnostics, never assume map is linear.
  reduced_gradient <- function(u) {
    if (!length(u)) return(numeric())
    if (is.null(grad)) return(.vpc_reduced_fd_gradient(reduced,u))
    par <- evaluate(u)$par
    if (is.null(par)) return(rep(NA_real_,length(u)))
    g <- tryCatch(grad(par),error=function(e) NULL)
    if (!is.numeric(g) || length(g)!=length(par) || any(!is.finite(g)))
      return(rep(NA_real_,length(u)))
    jac <- vapply(seq_along(u),function(i) {
      h <- 1e-5*max(1,abs(u[i])); a <- b <- u
      a[i] <- a[i]+h; b[i] <- b[i]-h
      pa <- tryCatch(problem$map(a),error=function(e) NULL)
      pb <- tryCatch(problem$map(b),error=function(e) NULL)
      if (!is.numeric(pa) || !is.numeric(pb) || length(pa)!=length(par) ||
          length(pb)!=length(par) || any(!is.finite(c(pa,pb))))
        return(rep(NA_real_,length(par)))
      (pa-pb)/(2*h)
    }, numeric(length(par)))
    as.numeric(crossprod(jac,g))
  }

  attempts <- list()
  for (i in seq_along(all_starts)) for (method in methods) {
    s <- all_starts[[i]]; tag <- paste0("start",i,"_",method)
    warnings <- character()
    t0 <- proc.time()[["elapsed"]]
    ans <- tryCatch(withCallingHandlers({
      if (!length(s)) list(par=s,convergence=0L,message="No free parameters")
      else if (method=="nlminb") {
        z <- stats::nlminb(s,reduced,control=list(iter.max=as.integer(maxit),
                                                    eval.max=as.integer(2*maxit)))
        list(par=z$par,convergence=z$convergence,message=z$message)
      } else {
        z <- stats::optim(s,reduced,method="BFGS",control=list(maxit=as.integer(maxit)))
        list(par=z$par,convergence=z$convergence,message=z$message)
      }
    },warning=function(w) { warnings <<- c(warnings,conditionMessage(w));
                            invokeRestart("muffleWarning") }),
    error=function(e) list(error=conditionMessage(e)))
    ans$elapsed_seconds <- proc.time()[["elapsed"]]-t0
    ans$warnings <- warnings
    if (is.null(ans$error)) {
      ev <- evaluate(ans$par)
      ans$objective <- ev$objective
      ans$violation <- ev$violation
      ans$valid <- is.finite(ev$objective)
      ans$optimizer_ok <- identical(as.integer(ans$convergence),0L)
      if (ans$valid) {
        g <- reduced_gradient(ans$par)
        ans$max_abs_reduced_gradient <- if (length(g)) {
          if (all(is.finite(g))) max(abs(g)) else Inf
        } else 0
        ans$stationary <- is.finite(ans$max_abs_reduced_gradient) &&
          ans$max_abs_reduced_gradient <= grad_tol
      } else { ans$stationary <- FALSE; ans$max_abs_reduced_gradient <- Inf }
      ans$numerically_acceptable <- ans$valid && ans$optimizer_ok && ans$stationary &&
        !is.na(ans$violation) && ans$violation<=feasibility_tol
    } else {
      ans$valid <- FALSE; ans$numerically_acceptable <- FALSE
    }
    attempts[[tag]] <- ans
  }
  valid_names <- names(attempts)[vapply(attempts,function(a) isTRUE(a$valid),logical(1))]
  if (!length(valid_names)) stop("No finite feasible constrained fit found.")
  acceptable_names <- valid_names[vapply(attempts[valid_names],function(a)
    isTRUE(a$numerically_acceptable),logical(1))]
  if (require_acceptable && !length(acceptable_names))
    stop("No numerically acceptable, independently validated constrained fit found.")
  # Keep both: best objective and best acceptable. A poorer acceptable
  # candidate does NOT establish the global restricted maximum.
  pick <- function(keys) keys[which.min(vapply(attempts[keys],`[[`,numeric(1),"objective"))]
  best_objective_name <- pick(valid_names)
  best_acceptable_name <- if (length(acceptable_names)) pick(acceptable_names) else NA_character_
  selected <- if (length(acceptable_names)) best_acceptable_name else best_objective_name
  best <- attempts[[selected]]
  full <- evaluate(best$par)$par
  obj_best <- attempts[[best_objective_name]]$objective
  obj_gap <- best$objective-obj_best
  diagnostics <- list(
    optimizer=selected, convergence_code=best$convergence,
    optimizer_ok=best$optimizer_ok,
    max_abs_reduced_gradient=best$max_abs_reduced_gradient,
    stationary=best$stationary,
    feasible=if (is.null(problem$validate)) NA else best$violation<=feasibility_tol,
    constraint_violation=best$violation,
    feasibility_tol=feasibility_tol,
    validation_provided=!is.null(problem$validate),
    numerically_acceptable=best$numerically_acceptable,
    n_attempts=length(attempts),
    best_objective_attempt=best_objective_name,
    best_objective=obj_best,
    best_acceptable_attempt=best_acceptable_name,
    selected_objective_gap=obj_gap,
    unresolved_better_candidate=obj_gap>1e-7*max(1,abs(obj_best)),
    global_optimality_verified=FALSE,
    boundary_complete=FALSE,
    covariance_strata_checked="not specified by generic parameter map"
  )
  structure(list(label=problem$label,par=full,free_par=best$par,
                 logLik=-best$objective,objective=best$objective,
                 diagnostics=diagnostics,attempts=attempts),
            class="glmmVPC_constrained_fit")
}

#' Compute observed likelihood-ratio statistic
#' @export
vpc_lrt_statistic <- function(full_fit,null_fit,negative_tol=1e-6) {
  ll1 <- full_fit$logLik; ll0 <- null_fit$logLik
  if (!is.numeric(ll1) || length(ll1)!=1L || !is.finite(ll1) ||
      !is.numeric(ll0) || length(ll0)!=1L || !is.finite(ll0))
    stop("Both fits must provide finite scalar logLik values.")
  if (!is.numeric(negative_tol) || length(negative_tol)!=1L ||
      !is.finite(negative_tol) || negative_tol<0)
    stop("negative_tol must be nonnegative and finite.")
  t <- 2*(ll1-ll0)
  if (t < -negative_tol)
    stop("Restricted log-likelihood exceeds full log-likelihood: check fits.")
  if (t < 0) warning("Small negative LRT due to numerical tolerance; returning zero.")
  max(0,t)
}

#' @export
print.glmmVPC_constrained_fit <- function(x,...) {
  cat("Constrained GLMM likelihood fit: ",x$label,"\n",sep="")
  cat("logLik: ",format(x$logLik,digits=9),
      " | acceptable: ",x$diagnostics$numerically_acceptable,
      " | global optimality verified: FALSE\n",sep="")
  invisible(x)
}
