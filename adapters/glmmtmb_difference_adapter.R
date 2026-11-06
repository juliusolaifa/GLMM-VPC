# Difference-contrast LRT adapter for the glmmTMB engine.
# The native glmmTMB theta parameterization covers Gamma > 0 only.
# This adapter also evaluates the vertex and available rank-one candidates,
# but does not certify the global restricted optimum across all PSD strata.
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
    sigma2 <- unname(family_params['sigma2'])
    if (length(sigma2) != 1L || !is.finite(sigma2) || sigma2 <= 0)
      stop('Gaussian residual variance sigma2 must be positive.')
    return(v/(v+sigma2))
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
                 c(sigma2=exp(disp))
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



#' Gaussian equality-null rank-one candidates
#'
#' For Gaussian VPC, equality at distinct x0,x1 is equivalent to
#' (u'z(x1))^2=(u'z(x0))^2. Thus u2=0 or
#' 2*u1+(x0+x1)*u2=0. Both lines are fitted exactly.
#' This reduction DOES NOT hold for non-Gaussian response-scale VPCs.
#' @export
fit_vpc_difference_rankone_gaussian <- function(
    fit, x0, x1, require_acceptable=TRUE) {
  if (!identical(fit$family,"gaussian"))
    stop("This exact two-direction reduction applies only to Gaussian VPC.")
  if (!requireNamespace("lme4",quietly=TRUE))
    stop("Install lme4.")
  if (!exists(".glmmvpc_glmmtmb_adapter",mode="function"))
    stop("Source fit_glmm.R first.")
  .vpc_glmmtmb_assert_fit(fit)
  if (x0==x1) stop("x0 and x1 must differ.")
  if (!is.data.frame(fit$data)) stop("Full fit must retain data.")
  bars <- lme4::findbars(fit$formula)
  if (length(bars)!=1L) stop("Expected one random intercept/slope term.")
  group <- fit$group
  if (!is.character(group) || length(group)!=1L ||
      !grepl("^[A-Za-z.][A-Za-z0-9._]*$",group))
    stop("Expected simple grouping column name.")
  # Reuse the actual random-effects design, including transformed slopes.
  re_formula <- stats::as.formula(call("~", bars[[1L]][[2L]]),
                                  env=environment(fit$formula))
  Xr <- stats::model.matrix(re_formula, data=fit$data)
  if (ncol(Xr)!=2L || colnames(Xr)[1L]!="(Intercept)")
    stop("Expected intercept and one random-slope design column.")
  slope <- Xr[,2L]
  if (any(!is.finite(slope))) stop("Random slope design is not finite.")
  directions <- list(intercept_only=c(1,0),
                     centered_slope=c(-(x0+x1)/2,1))
  candidates <- lapply(names(directions),function(nm) {
    tryCatch({
    u <- directions[[nm]]
    dat <- fit$data
    col <- ".glmmvpc_difference_re"
    if (col %in% names(dat)) stop("Reserved column already exists.")
    dat[[col]] <- u[1]+u[2]*slope
    fixed <- lme4::nobars(fit$formula)
    re_term <- call("|",call("+",0,as.name(col)),as.name(group))
    f <- stats::as.formula(call("~",fixed[[2]],
       call("+",fixed[[3]],call("(",re_term))),
       env=environment(fit$formula))
    ans <- .glmmvpc_glmmtmb_adapter(f,dat,fit$family,fit$control)
    .vpc_glmmtmb_compare_rows(fit$raw_fit,ans$fit)
    ll <- as.numeric(ans$diagnostics$logLik)
    check <- .vpc_glmmtmb_objective_check(
      ans$fit,ll,.vpc_glmmtmb_ll_tol(fit))
    if (!check$consistent) stop("Objective mismatch.")
    vc <- as.matrix(glmmTMB::VarCorr(ans$fit)$cond[[group]])
    if (!identical(dim(vc),c(1L,1L)) || !is.finite(vc[1,1]) ||
        vc[1,1]<0) stop("Invalid scalar random-effect variance.")
    G <- unname(vc[1,1])*tcrossprod(u)
    beta <- glmmTMB::fixef(ans$fit)$cond
    fp <- .glmmvpc_extract_null_family_params(ans$fit,fit$family)
    d <- vpc_response_difference(beta,G,fit$family,fp,x0,x1)
    if (!is.finite(d) || abs(d)>1e-8)
      stop("Gaussian rank-one equality check failed.")
    dg <- ans$diagnostics
    dg$feasible <- TRUE
    dg$constraint_violation <- abs(d)
    dg$rankone_direction <- u
    dg$selected_stratum <- paste0("rankone_",nm)
    dg$boundary_complete <- FALSE
    dg$global_optimality_verified <- FALSE
    structure(list(label="difference Gaussian rank-one equality",
      alternative="two.sided",x0=x0,x1=x1,contrast=d,
      beta=beta,Gamma=G,family_params=fp,
      logLik=ll,objective=check$objective,engine="glmmTMB",
      raw_fit=ans$fit,formula=f,diagnostics=dg,
      attempts=ans$attempts),class="glmmVPC_constrained_fit")
    }, error=function(e) e)
  })
  names(candidates) <- names(directions)
  candidates
}


#' Enumerate rank-one equality directions at fixed nuisance parameters
#'
#' A numerical geometry helper, not a likelihood optimizer. The roots
#' may be incomplete when the contrast touches zero without changing sign.
#' @export
vpc_difference_rankone_roots <- function(
    beta, tau2, family, family_params, x0, x1,
    fixed_design=NULL, grid_size=720L, tol=1e-9) {
  if (!is.numeric(tau2) || length(tau2)!=1L ||
      !is.finite(tau2) || tau2<=0)
    stop("tau2 must be strictly positive.")
  if (!is.numeric(grid_size) || length(grid_size)!=1L ||
      grid_size<16L) stop("grid_size must be at least 16.")
  if (x0==x1) stop("x0 and x1 must differ.")
  grid <- seq(0,pi,length.out=as.integer(grid_size)+1L)
  difference <- function(angle) {
    u <- c(cos(angle),sin(angle))
    vpc_response_difference(
      beta,tau2*tcrossprod(u),family,family_params,
      x0,x1,fixed_design)
  }
  vals <- vapply(grid,difference,numeric(1))
  if (any(!is.finite(vals))) stop("Nonfinite contrast along rank-one circle.")
  roots <- numeric(0)
  for (i in seq_len(length(grid)-1L)) {
    if (abs(vals[i])<=tol) roots <- c(roots,grid[i])
    if (vals[i]*vals[i+1L]<0) {
      roots <- c(roots,stats::uniroot(
        difference,grid[c(i,i+1L)],tol=tol)$root)
    }
  }
  if (abs(vals[length(vals)])<=tol) roots <- c(roots,pi)
  roots <- sort(roots %% pi)
  roots <- roots[!duplicated(round(roots,10))]
  list(angles=roots,
       directions=lapply(roots,function(a)c(cos(a),sin(a))),
       contrast=vapply(roots,difference,numeric(1)),
       grid_size=grid_size,
       exhaustive=FALSE,
       warning=paste(
         "Sign-change enumeration can miss tangent or multiple roots;",
         "nuisance parameters and variance were held fixed;",
         "this is not a restricted maximum likelihood fit."))
}

#' Fit available strata of the response-scale VPC difference null
#'
#' Numerical candidates are compared at the vertex, interior, and available
#' rank-one strata. The rank-one nonlinear search is experimental and the
#' global restricted optimum is NOT certified.
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
    d$boundary_complete <- FALSE
    d$global_optimality_verified <- FALSE
    d$warning <- "LRT is zero conditional on unrestricted fit being the global MLE."
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
  rankone <- list()
  if (identical(fit$family, "gaussian") &&
      identical(alternative, "two.sided")) {
    rankone <- tryCatch(
      fit_vpc_difference_rankone_gaussian(
        fit, x0, x1, require_acceptable=require_acceptable),
      error=function(e) list(rankone_error=e))
  }
  # Non-Gaussian equality requires nonlinear rank-one fitting. This
  # experimental search is not a global optimizer; failures are retained.
  if (!identical(fit$family, "gaussian")) {
    rankone <- list(rankone_nonlinear = tryCatch(
      .vpc_difference_rankone_nonlinear(
        fit, x0, x1, alternative = alternative,
        fixed_design = fixed_design, maxeval = min(as.integer(maxeval), 400L),
        feasibility_tol = feasibility_tol,
        require_acceptable = require_acceptable),
      error = function(e) e))
  }
  candidates <- c(list(interior=interior,vertex=vertex),rankone)
  ok <- vapply(candidates, valid, logical(1))
  if (!any(ok)) stop("No usable restricted covariance stratum.")
  lls <- vapply(names(candidates), function(k)
    if (ok[[k]]) candidates[[k]]$logLik else NA_real_, numeric(1))
  chosen <- names(which.max(replace(lls, is.na(lls), -Inf)))
  out <- candidates[[chosen]]
  out$diagnostics$selected_stratum <- chosen
  out$diagnostics$stratum_logLik <- lls
  out$diagnostics$covariance_strata_checked <- names(ok)[ok]
  out$diagnostics$stratum_errors <- lapply(candidates, function(z)
    if (inherits(z, "error")) conditionMessage(z) else NULL)
  # Coverage of covariance strata is not a certificate of a global optimum.
  out$diagnostics$boundary_complete <- FALSE
  out$diagnostics$strata_enumerated <-
    identical(fit$family,"gaussian") &&
    identical(alternative,"two.sided") &&
    all(ok[c("interior","vertex","intercept_only","centered_slope")])
  out$diagnostics$global_optimality_verified <- FALSE
  out$diagnostics$warning <-
    "Numerical candidates only; global restricted maximum and complete PSD alternative not certified."
  out
}

#' Observed difference LRT from provisional stratified restricted fitting
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
       calibration = "not computed; global restricted optimum not certified")
}
\n\n# ---- Internal experimental rank-one nonlinear likelihood candidate ----\n# Experimental rank-one constrained marginal likelihood optimization.
# Source after fit_glmm.R and glmmtmb_difference_adapter.R.
# Uses glmmTMB's Laplace approximation with ALL outer parameters fixed
# at each evaluation, allowing the random-effects integral to be recomputed.
# Extremely expensive; not a production-certified global optimizer.

.vpc_rankone_family <- function(family) switch(family,
  gaussian=stats::gaussian(), poisson=stats::poisson(),
  nbinom=glmmTMB::nbinom2(link="log"),
  cpg=glmmTMB::tweedie(link="log"),
  stop("Unsupported family"))

.vpc_rankone_fixed_model <- function(fit, angle, par, fixed_design=NULL) {
  if (!requireNamespace("glmmTMB",quietly=TRUE) ||
      !requireNamespace("lme4",quietly=TRUE))
    stop("Install glmmTMB and lme4.")
  if (length(lme4::findbars(fit$formula))!=1L)
    stop("Requires one random intercept-and-slope term.")
  bar <- lme4::findbars(fit$formula)[[1L]]
  re <- stats::model.matrix(stats::as.formula(call("~",bar[[2L]]),
    env=environment(fit$formula)),fit$data)
  if (ncol(re)!=2L || colnames(re)[1L]!="(Intercept)")
    stop("Expected two-column random design, intercept first.")
  u <- c(cos(angle),sin(angle))
  dat <- fit$data
  nm <- ".glmmvpc_rankone_design"
  if (nm %in% names(dat)) stop("Reserved design column collision.")
  dat[[nm]] <- drop(re %*% u)
  grp <- fit$group
  if (!is.character(grp) || length(grp)!=1L ||
      !grepl("^[A-Za-z.][A-Za-z0-9._]*$",grp))
    stop("Expected simple grouping variable.")
  fixed <- lme4::nobars(fit$formula)
  term <- call("(",call("|",call("+",0L,as.name(nm)),as.name(grp)))
  formula <- stats::as.formula(call("~",fixed[[2L]],
    call("+",fixed[[3L]],term)),env=environment(fit$formula))
  # Outer vector is (beta, log_sd, log_sigma2/size/phi, power_logit).
  nb <- length(fit$beta)
  if (length(par)!=nb+1L+length(fit$family_params))
    stop("Outer parameter length mismatch.")
  beta <- par[seq_len(nb)]
  theta <- par[nb+1L]
  disp <- if (fit$family %in% c("gaussian","nbinom","cpg"))
    par[nb+2L] else numeric()
  psi <- if (fit$family=="cpg") par[nb+3L] else numeric()
  start <- list(beta=unname(beta),theta=theta)
  map <- list(beta=factor(rep(NA_integer_,nb)),theta=factor(NA_integer_))
  if (length(disp)) {
    start$betadisp <- disp
    map$betadisp <- factor(NA_integer_)
  }
  if (length(psi)) {
    start$psi <- psi
    map$psi <- factor(NA_integer_)
  }
  raw <- suppressWarnings(glmmTMB::glmmTMB(
    formula,data=dat,family=.vpc_rankone_family(fit$family),
    start=start,map=map))
  ll <- as.numeric(stats::logLik(raw))
  if (!is.finite(ll)) stop("Nonfinite rank-one marginal likelihood.")
  fp <- switch(fit$family,
    gaussian=c(sigma2=exp(disp)),poisson=numeric(),
    nbinom=c(theta=exp(disp)),
    cpg=c(phi=exp(disp),power=1+stats::plogis(psi)))
  list(logLik=ll,beta=stats::setNames(beta,names(fit$beta)),
       Gamma=exp(2*theta)*tcrossprod(u),family_params=fp,
       raw_fit=raw,angle=angle,formula=formula)
}

#' Experimental constrained rank-one GLMM difference fit
#'
#' This is a computational prototype. Rebuilds the Laplace likelihood at
#' each trial angle and parameter vector. Does not certify a global optimum.
#' @export
.vpc_difference_rankone_nonlinear <- function(
    fit,x0,x1,alternative=c("two.sided","greater","less"),
    fixed_design=NULL,maxeval=120L,feasibility_tol=1e-6,
    require_acceptable=TRUE) {
  alternative <- match.arg(alternative)
  if (!inherits(fit,"glmmVPC_fit") || fit$engine!="glmmTMB")
    stop("Requires glmmTMB glmmVPC_fit.")
  if (x0==x1) stop("Distinct covariates required.")
  if (!requireNamespace("nloptr",quietly=TRUE)) stop("Install nloptr.")
  if (require_acceptable && !isTRUE(fit$diagnostics$numerically_acceptable))
    stop("Full fit is not numerically acceptable.")
  nb <- length(fit$beta)
  eig <- eigen(fit$Gamma,symmetric=TRUE)
  u0 <- eig$vectors[,1L]
  angle0 <- atan2(u0[2L],u0[1L]) %% pi
  sd0 <- sqrt(max(eig$values[1L],1e-5))
  fp <- switch(fit$family,gaussian=log(fit$family_params[["sigma2"]]),
    poisson=numeric(),nbinom=log(fit$family_params[["theta"]]),
    cpg=c(log(fit$family_params[["phi"]]),
          stats::qlogis(fit$family_params[["power"]]-1)))
  start <- c(unname(fit$beta),log(sd0),fp,angle0)
  # The angle is periodic with period pi. We leave it unconstrained
  # and use multiple starts to avoid an artificial seam at 0/pi.
  cache <- new.env(parent=emptyenv())
  evaluate <- function(x) {
    key <- paste(formatC(x,digits=14,format="fg"),collapse="|")
    if (exists(key,envir=cache,inherits=FALSE))
      return(get(key,envir=cache))
    result <- tryCatch({
      z <- .vpc_rankone_fixed_model(fit,tail(x,1L),head(x,-1L),fixed_design)
      d <- vpc_response_difference(z$beta,z$Gamma,fit$family,
        z$family_params,x0,x1,fixed_design)
      list(ok=TRUE,fit=z,difference=d,nll=-z$logLik)
    },error=function(e) list(ok=FALSE,error=conditionMessage(e),
                             difference=NA_real_,nll=1e30))
    assign(key,result,envir=cache)
    result
  }
  objective <- function(x) evaluate(x)$nll
  constraint <- function(x) {
    z <- evaluate(x)
    if (!z$ok) return(1e6)
    switch(alternative,two.sided=z$difference,
      greater=z$difference,less=-z$difference)
  }
  # COBYLA supports equality and inequality constraints without numerical
  # gradients, but can be slow and may fail for tight feasibility tolerances.
  starts <- lapply(c(0,pi/4,pi/2,3*pi/4),function(a) {
    x <- start; x[length(x)] <- a; x
  })
  attempts <- lapply(starts,function(x) tryCatch(nloptr::nloptr(
    x0=x,eval_f=objective,
    eval_g_eq=if(alternative=="two.sided") constraint else NULL,
    eval_g_ineq=if(alternative!="two.sided") constraint else NULL,
    opts=list(algorithm="NLOPT_LN_COBYLA",maxeval=as.integer(maxeval),
              xtol_rel=1e-6)),error=function(e)e))
  valid <- vapply(attempts,function(a) {
    if (inherits(a,"error") || is.null(a$solution)) return(FALSE)
    z <- evaluate(a$solution)
    z$ok && is.finite(z$difference) &&
      (if(alternative=="two.sided") abs(z$difference)<=feasibility_tol
       else if(alternative=="greater") z$difference<=feasibility_tol
       else z$difference>=-feasibility_tol)
  },logical(1))
  if (!any(valid)) stop("No feasible nonlinear rank-one fit; attempts=",
                        length(attempts),". This is NOT proof of infeasibility.")
  ii <- which(valid)
  best <- ii[which.min(vapply(attempts[ii],function(a)
    evaluate(a$solution)$nll,numeric(1)))]
  opt <- attempts[[best]]
  z <- evaluate(opt$solution)
  diagnostics <- list(feasible=TRUE,
    numerically_acceptable=isTRUE(opt$status>0),
    constraint_violation=if(alternative=="two.sided") abs(z$difference)
      else if(alternative=="greater") max(0,z$difference)
      else max(0,-z$difference),
    optimizer="COBYLA / repeated glmmTMB Laplace evaluations",
    optimizer_status=opt$status,selected_stratum="rankone_nonlinear",
    boundary_complete=FALSE,global_optimality_verified=FALSE,
    warning="Experimental nested likelihood; verify fixed-parameter map and convergence.")
  structure(c(list(label="nonlinear rank-one difference",
    alternative=alternative,x0=x0,x1=x1,contrast=z$difference,
    logLik=z$fit$logLik,objective=z$nll,engine="glmmTMB",
    diagnostics=diagnostics,attempts=attempts),
    z$fit[c("beta","Gamma","family_params","raw_fit","angle","formula")]),
    class="glmmVPC_constrained_fit")
}
