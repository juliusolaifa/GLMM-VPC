# Observed omnibus likelihood-ratio statistic for H0: Gamma = 0.
# The restricted GLMM is fitted by fit_vpc_omnibus_null() in
# glmmtmb_likelihood_adapter.R. Calibration is owned by lrt_calibration.R.

#' Compute the observed omnibus LRT
#' @param fit Fitted glmmVPC_fit object from the glmmTMB engine.
#' @param require_acceptable Require acceptable numerical diagnostics for both fits.
#' @param negative_tol Tolerance for small negative LRT values from optimization.
#' @return Uncalibrated statistic, full/null fits, log-likelihoods, diagnostics.
#' @export
vpc_omnibus_lrt <- function(fit, require_acceptable = TRUE,
                            negative_tol = 1e-6) {
  if (!inherits(fit, "glmmVPC_fit") || !identical(fit$engine, "glmmTMB"))
    stop("fit must be a glmmTMB glmmVPC_fit object.", call. = FALSE)
  if (!is.logical(require_acceptable) || length(require_acceptable) != 1L ||
      is.na(require_acceptable))
    stop("require_acceptable must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(negative_tol) || length(negative_tol) != 1L ||
      !is.finite(negative_tol) || negative_tol < 0)
    stop("negative_tol must be a finite nonnegative scalar.", call. = FALSE)
  if (require_acceptable && !isTRUE(fit$diagnostics$numerically_acceptable))
    stop("Full fit is not numerically acceptable; inspect fit$diagnostics.",
         call. = FALSE)

  null_fit <- fit_vpc_omnibus_null(fit,
                                   require_acceptable = require_acceptable)
  stat <- vpc_lrt_statistic(fit, null_fit, negative_tol = negative_tol)
  list(statistic = stat, full_fit = fit, null_fit = null_fit,
       full_logLik = fit$logLik, null_logLik = null_fit$logLik,
       calibration = "not computed",
       diagnostics = list(full = fit$diagnostics,
                          null = null_fit$diagnostics))
}
