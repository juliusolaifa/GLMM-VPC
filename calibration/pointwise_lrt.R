# Pointwise zero-VPC likelihood-ratio calibration for a 2x2 PSD covariance.
# Depends on cone_geometry.R, cone_metric.R, cone_projection.R.
# All information matrices are efficient information for gamma=(g11,g12,g22).

#' Pointwise zero-variance null geometry
#' @param x0 Finite covariate value.
#' @param information Positive definite efficient 3x3 information matrix.
#' @return List containing covariance and whitened null-ray directions.
#' @export
pointwise_null_geometry <- function(x0, information) {
  x0 <- .cg_scalar(x0, "x0")
  I <- .cm_information(information)
  list(x0 = x0, gamma_ray = pointwise_null_ray(x0, "gamma"),
       whitened_ray = whitened_pointwise_ray(x0, I), information = I)
}

#' Gaussian-limit pointwise LRT statistic
#'
#' This computes distance(y,C0)^2 - distance(y,C)^2 where C is the
#' whitened PSD cone and C0 is the whitened pointwise null ray.
#' @param y Whitened Gaussian-limit observation, length three.
#' @param x0 Covariate value defining the null.
#' @param information Efficient information matrix.
#' @param ... Numerical controls passed to project_whitened_psd().
#' @return A nonnegative scalar.
#' @export
pointwise_lrt_limit_statistic <- function(y, x0, information, ...) {
  pointwise_projection_statistic(y, x0, information, ...)
}

#' Calibrate a pointwise zero-VPC likelihood-ratio statistic
#'
#' Uses the maximum of the vertex and nonzero-ray endpoint
#' survival probabilities and critical values.
#'
#' The vertex distribution is evaluated by deterministic
#' angular quadrature. The nonzero-ray endpoint follows
#' 0.5*chisq(1) + 0.5*chisq(2).
#'
#' @param statistic Nonnegative observed LRT statistic.
#' @param x0 Covariate value defining the pointwise null.
#' @param information Positive-definite efficient information
#'   matrix for gamma=(g11,g12,g22), evaluated at the
#'   pointwise restricted fit.
#' @param alpha Significance level in (0,1).
#' @param ... Numerical quadrature controls passed to
#'   calibrate_pointwise_angular().
#' @return List with unified and component calibration results.
#' @export
calibrate_pointwise_lrt <- function(
    statistic,
    x0,
    information,
    alpha = 0.05,
    ...
) {

  statistic <- .cg_scalar(statistic, "statistic")

  if (statistic < 0) {
    stop(
      "statistic must be nonnegative.",
      call. = FALSE
    )
  }

  x0 <- .cg_scalar(x0, "x0")
  information <- .cm_information(information)

  if (
    !is.numeric(alpha) ||
    length(alpha) != 1L ||
    !is.finite(alpha) ||
    alpha <= 0 ||
    alpha >= 1
  ) {
    stop(
      "alpha must lie strictly between 0 and 1.",
      call. = FALSE
    )
  }

  # Unified deterministic endpoint calibration
  cal <- calibrate_pointwise_angular(
    statistic = statistic,
    x0 = x0,
    information = information,
    alpha = alpha,
    ...
  )

  # Required output fields
  required <- c(
    "p_value",
    "critical_value",
    "reject",
    "vertex_p_value",
    "ray_p_value",
    "vertex_critical",
    "ray_critical"
  )

  missing <- setdiff(required, names(cal))

  if (length(missing)) {
    stop(
      paste(
        "Angular calibration did not return:",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Check consistency with the unified envelope
  expected_p <- max(
    cal$vertex_p_value,
    cal$ray_p_value
  )

  expected_critical <- max(
    cal$vertex_critical,
    cal$ray_critical
  )

  if (
    !isTRUE(all.equal(
      cal$p_value,
      expected_p,
      tolerance = 1e-8
    )) ||
    !isTRUE(all.equal(
      cal$critical_value,
      expected_critical,
      tolerance = 1e-8
    ))
  ) {
    stop(
      "Inconsistent pointwise endpoint calibration.",
      call. = FALSE
    )
  }

  if (!isTRUE(all.equal(cal$reject, statistic > cal$critical_value))) {
    stop("Inconsistent pointwise rejection decision.", call. = FALSE)
  }

  list(
    statistic = statistic,
    x0 = x0,
    alpha = alpha,
    information = information,

    p_value = cal$p_value,
    critical_value = cal$critical_value,
    reject = cal$reject,

    vertex_p_value = cal$vertex_p_value,
    ray_p_value = cal$ray_p_value,

    vertex_critical = cal$vertex_critical,
    ray_critical = cal$ray_critical,

    method = "deterministic_angular_endpoint_envelope"
  )
}


#' Observed pointwise zero-VPC LRT (closed null ray, including vertex)
#'
#' The null fitter in glmmtmb_pointwise_adapter.R compares the rank-one
#' ray and Gamma=0. This routine does not compute a p-value: the angular
#' reference is supplied separately by calibrate_pointwise_lrt().
#' @param fit Fitted glmmVPC_fit from the glmmTMB engine.
#' @param x0 Finite pointwise covariate value.
#' @param require_acceptable Require acceptable full and restricted fits.
#' @param negative_tol Allowed numerical likelihood-ordering discrepancy.
#' @return Uncalibrated LRT result and restricted-fit diagnostics.
#' @export
vpc_pointwise_lrt <- function(fit, x0, require_acceptable = TRUE,
                              negative_tol = 1e-6) {
  if (!inherits(fit, "glmmVPC_fit") || !identical(fit$engine, "glmmTMB"))
    stop("fit must be a glmmTMB glmmVPC_fit object.", call. = FALSE)
  if (!is.numeric(x0) || length(x0) != 1L || !is.finite(x0))
    stop("x0 must be a finite numeric scalar.", call. = FALSE)
  if (!is.logical(require_acceptable) || length(require_acceptable) != 1L ||
      is.na(require_acceptable))
    stop("require_acceptable must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(negative_tol) || length(negative_tol) != 1L ||
      !is.finite(negative_tol) || negative_tol < 0)
    stop("negative_tol must be finite and nonnegative.", call. = FALSE)
  if (require_acceptable && !isTRUE(fit$diagnostics$numerically_acceptable))
    stop("Full fit is not numerically acceptable.", call. = FALSE)

  null_fit <- fit_vpc_pointwise_null(fit, x0,
                                     require_acceptable = require_acceptable)
  if (require_acceptable && !isTRUE(null_fit$diagnostics$boundary_complete))
    stop("Pointwise null fit did not cover both vertex and ray strata.",
         call. = FALSE)
  stat <- vpc_lrt_statistic(fit, null_fit, negative_tol = negative_tol)
  list(statistic = stat, x0 = x0, full_fit = fit, null_fit = null_fit,
       full_logLik = fit$logLik, null_logLik = null_fit$logLik,
       calibration = "not computed",
       diagnostics = list(full = fit$diagnostics,
                          null = null_fit$diagnostics))
}
