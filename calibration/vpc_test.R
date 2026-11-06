# Unified public likelihood-ratio interface for glmmVPC.
# Source the appropriate vpc_*_lrt adapters separately.
# Numerical calibrations are registered as providers, not fabricated here.
#
# control$calibrators is a named list with entries:
#   omnibus$verified_uniform / $conditional_uniform / $gaussian_local / $pointwise
#   pointwise$verified_uniform / $conditional_uniform / $gaussian_local / $pointwise
#   difference$two.sided$... / $greater$... / $less$...
# For difference, the intended gaussian_local providers are the improved
# U_Q envelope (two.sided) and directional cone-distance bound (greater/less).
#
# Each provider must be a function(fit, null_fit, statistic, alpha,
# information, x0, x1, alternative, control), returning a list with
# p.value, critical.value, method, validity, assumptions (character vector).
# validity must match its registry key. This declaration is not a proof.

.vpc_test_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop(name, " must be one finite numeric value.")
  invisible(x)
}

.vpc_test_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE.")
  invisible(x)
}

.vpc_test_empty_calibration <- function(reason = "No calibration provider connected") {
  list(p.value = NA_real_, critical.value = NA_real_, reject = NA,
       method = "not calibrated", validity = "not_established",
       assumptions = character(), reason = reason)
}

.vpc_test_provider_group <- function(control, test, alternative) {
  registry <- control$calibrators
  if (is.null(registry)) return(list())
  if (!is.list(registry)) stop("control$calibrators must be a named list.")
  group <- registry[[test]]
  if (is.null(group)) return(list())
  if (!is.list(group)) stop("control$calibrators[[test]] must be a list.")
  if (test == "difference") {
    group <- group[[alternative]]
    if (is.null(group)) return(list())
    if (!is.list(group))
      stop("Difference calibrators must be grouped by alternative.")
  }
  group
}

.vpc_test_calibrate <- function(test, requested, fit, observed, alpha,
                                information, x0, x1, alternative, control) {
  if (requested == "none")
    return(.vpc_test_empty_calibration("Calibration explicitly disabled"))
  group <- .vpc_test_provider_group(control, test, alternative)
  order <- switch(requested,
    auto = c("verified_uniform", "conditional_uniform", "gaussian_local", "pointwise"),
    uniform = "verified_uniform",
    pointwise = "pointwise")
  selected <- NULL
  for (key in order) {
    if (is.function(group[[key]])) { selected <- key; break }
  }
  if (is.null(selected)) {
    if (requested != "auto")
      stop("No '", requested, "' provider connected for ", test,
           if (test == "difference") paste0(" (", alternative, ")") else "", ".")
    return(.vpc_test_empty_calibration())
  }
  # Only a separately justified and explicitly enabled conditional theorem
  # can be described as applicable to the user's GLMM setting.
  if (selected == "conditional_uniform" &&
      !isTRUE(control$conditional_assumptions_verified)) {
    if (requested == "auto") {
      selected <- NULL
      for (key in c("gaussian_local", "pointwise")) {
        if (is.function(group[[key]])) { selected <- key; break }
      }
      if (is.null(selected))
        return(.vpc_test_empty_calibration(
          "Conditional-uniform assumptions not verified; no fallback connected"))
    }
  }
  provider <- group[[selected]]
  ans <- provider(fit = fit, null_fit = observed$null_fit,
                  statistic = observed$statistic, alpha = alpha,
                  information = information, x0 = x0, x1 = x1,
                  alternative = alternative, control = control)
  if (!is.list(ans)) stop("Calibration provider must return a list.")
  if (selected == "conditional_uniform" &&
      !isTRUE(control$conditional_assumptions_verified))
    stop("Conditional-uniform calibration requires verified assumptions.")
  if (!identical(ans$validity, selected))
    stop("Calibration provider validity must equal its registry key: '", selected, "'.")
  if (!is.character(ans$method) || length(ans$method) != 1L ||
      is.na(ans$method) || !nzchar(ans$method))
    stop("Calibration provider must return a nonempty method string.")
  if (!is.character(ans$assumptions) || anyNA(ans$assumptions))
    stop("Calibration provider must return assumptions as a character vector.")
  .vpc_test_scalar(ans$p.value, "Calibration p.value")
  .vpc_test_scalar(ans$critical.value, "Calibration critical.value")
  if (ans$p.value < 0 || ans$p.value > 1 || ans$critical.value < 0)
    stop("Invalid calibration p-value or critical value.")
  # Strict rejection respects atoms at zero and conservative upper tails.
  ans$reject <- observed$statistic > ans$critical.value
  ans$requested <- requested
  ans$assumptions_verified <- if (selected == "conditional_uniform")
    isTRUE(control$conditional_assumptions_verified) else NA
  ans
}

#' Unified VPC likelihood-ratio testing
#'
#' @param fit A glmmVPC_fit object returned by fit_glmm().
#' @param test One of "omnibus", "pointwise", "difference".
#' @param x0 Covariate value for pointwise/difference tests.
#' @param x1 Second covariate value for difference tests.
#' @param alternative "two.sided", "greater", "less" for difference.
#' @param alpha Nominal size in (0,1).
#' @param calibration "auto" (best available eligible provider), "uniform"
#'   (requires verified_uniform), "pointwise" (explicit fixed-parameter
#'   calibration), or "none". Gaussian-local conservatism is not a proof
#'   of GLMM-uniform size control.
#' @param information Optional information object passed to provider. Providers
#'   must check parameter ordering, evaluation point, and whether efficient
#'   covariance information J or full information I is required.
#' @param control Named list. Supported entries: require_acceptable (TRUE),
#'   fixed_design, calibrators (provider registry),
#'   conditional_assumptions_verified (FALSE). Other entries pass to providers.
#' @return glmmVPC_test with standardized calibration metadata.
#' @export
vpc_test <- function(fit,
                     test = c("omnibus", "pointwise", "difference"),
                     x0 = NULL, x1 = NULL,
                     alternative = "two.sided",
                     alpha = 0.05,
                     calibration = c("auto", "uniform", "pointwise", "none"),
                     information = NULL,
                     control = list()) {
  test <- match.arg(test)
  calibration <- match.arg(calibration)
  if (!inherits(fit, "glmmVPC_fit")) stop("fit must be a glmmVPC_fit.")
  if (!is.list(control)) stop("control must be a list.")
  .vpc_test_scalar(alpha, "alpha")
  if (alpha <= 0 || alpha >= 1) stop("alpha must lie in (0,1).")
  acceptable <- control$require_acceptable
  if (is.null(acceptable)) acceptable <- TRUE
  .vpc_test_flag(acceptable, "control$require_acceptable")
  if (!is.null(control$conditional_assumptions_verified))
    .vpc_test_flag(control$conditional_assumptions_verified,
                   "control$conditional_assumptions_verified")

  if (test == "omnibus") {
    if (!is.null(x0) || !is.null(x1) || !identical(alternative, "two.sided"))
      stop("Omnibus accepts neither covariate values nor directional alternatives.")
    observed <- vpc_omnibus_lrt(fit, require_acceptable = acceptable)
  } else if (test == "pointwise") {
    .vpc_test_scalar(x0, "x0")
    if (!is.null(x1) || !identical(alternative, "two.sided"))
      stop("Pointwise requires x0 only; alternative must be 'two.sided'.")
    observed <- vpc_pointwise_lrt(fit, x0, require_acceptable = acceptable)
  } else {
    .vpc_test_scalar(x0, "x0")
    .vpc_test_scalar(x1, "x1")
    if (x0 == x1) stop("x0 and x1 must differ.")
    alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
    # The older difference_lrt.R defines vpc_difference_lrt(full_fit,
    # restricted, negative_tol), which is NOT a fitting/dispatch wrapper.
    # Fail explicitly rather than silently passing x0 as a restricted fit.
    if (!exists("vpc_difference_lrt", mode = "function") ||
        !all(c("x0", "x1", "alternative", "fixed_design",
               "require_acceptable") %in% names(formals(vpc_difference_lrt))))
      stop("Difference test dispatcher is not connected: the current ",
           "vpc_difference_lrt() accepts a pre-fitted restricted model, ",
           "not (fit, x0, x1). Integrate the canonical difference adapter first.")
    observed <- vpc_difference_lrt(
      fit, x0, x1, alternative = alternative,
      fixed_design = control$fixed_design,
      require_acceptable = acceptable)
    complete <- isTRUE(observed$diagnostics$boundary_complete) ||
      isTRUE(observed$null_fit$diagnostics$boundary_complete)
    if (!complete && calibration != "none")
      stop("Difference null fit is not certified boundary-complete. ",
           "Calibrated inference is blocked; use calibration='none' ",
           "for exploratory statistics until the full null optimizer is connected.")
  }
  .vpc_test_scalar(observed$statistic, "LRT statistic")
  if (observed$statistic < 0)
    stop("Negative LRT statistic: check nested fit objectives and convergence; do not silently truncate it.")
  cal <- .vpc_test_calibrate(test, calibration, fit, observed, alpha,
                             information, x0, x1, alternative, control)
  structure(list(
    test = test, x0 = x0, x1 = x1,
    alternative = if (test == "difference") alternative else NULL,
    statistic = observed$statistic,
    estimated_contrast = observed$estimated_contrast,
    p.value = cal$p.value, critical.value = cal$critical.value,
    reject = cal$reject, alpha = alpha,
    calibration_requested = calibration,
    calibration = list(method = cal$method, validity = cal$validity,
                       assumptions = cal$assumptions,
                       details = cal),
    # Compatibility convenience fields for existing plotting code.
    calibration_validity = cal$validity,
    information = information,
    full_fit = if (is.null(observed$full_fit)) fit else observed$full_fit,
    null_fit = observed$null_fit,
    diagnostics = observed$diagnostics
  ), class = "glmmVPC_test")
}

#' @export
print.glmmVPC_test <- function(x, ...) {
  cat("glmmVPC test: ", x$test, "\n", sep = "")
  if (!is.null(x$alternative)) cat("Alternative: ", x$alternative, "\n", sep = "")
  cat("LRT = ", format(x$statistic, digits = 8), "\n", sep = "")
  if (is.na(x$p.value)) {
    reason <- x$calibration$details$reason
    if (is.null(reason)) reason <- "No valid calibration available"
    cat("Calibration unavailable: ", reason, "\n", sep = "")
  } else {
    cat("p-value = ", format(x$p.value, digits = 8),
        "; critical = ", format(x$critical.value, digits = 8),
        "; reject = ", x$reject, "\n", sep = "")
    cat("Method: ", x$calibration$method, "\n", sep = "")
    cat("Validity: ", x$calibration$validity, "\n", sep = "")
  }
  invisible(x)
}
