# Chi-bar-square weights for testing H0: Gamma = 0,
# where Gamma is a 2x2 positive-semidefinite
# random-effects covariance matrix.
#
# Requires a 3x3 positive-definite efficient
# Fisher information matrix for
# gamma = (gamma11, gamma12, gamma22).


.validate_chibar_information <- function(information) {
  if (!is.matrix(information) || !is.numeric(information) ||
      !identical(dim(information), c(3L, 3L)) ||
      any(!is.finite(information))) {
    stop("information must be a finite numeric 3x3 matrix.", call. = FALSE)
  }
  if (!isTRUE(all.equal(information, t(information), tolerance = 1e-10))) {
    stop("information must be symmetric.", call. = FALSE)
  }
  information <- (information + t(information)) / 2
  if (inherits(try(chol(information), silent = TRUE), "try-error")) {
    stop("information must be positive definite.", call. = FALSE)
  }
  invisible(information)
}

.validate_chibar_quadrature <- function(rel.tol, subdivisions) {
  if (!is.numeric(rel.tol) || length(rel.tol) != 1L ||
      !is.finite(rel.tol) || rel.tol <= 0 || rel.tol >= 1) {
    stop("rel.tol must be a finite number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(subdivisions) || length(subdivisions) != 1L ||
      !is.finite(subdivisions) || subdivisions < 1 ||
      subdivisions != floor(subdivisions) ||
      subdivisions > .Machine$integer.max) {
    stop("subdivisions must be a positive integer.", call. = FALSE)
  }
  invisible(TRUE)
}


.matrix_sqrt_spd <- function(M, inverse = FALSE) {
  E <- base::eigen((M + t(M)) / 2, symmetric = TRUE)
  if (any(E$values <= 0)) stop("M must be positive definite.", call. = FALSE)
  d <- if (inverse) 1 / sqrt(E$values) else sqrt(E$values)
  base::sweep(E$vectors, 2L, d, `*`) %*% t(E$vectors)
}

.assemble_chibar_result <- function(w0, w3, method,
                                   err0 = NA_real_, err3 = NA_real_,
                                   se0 = NA_real_, se3 = NA_real_,
                                   nsim = NA_integer_) {
  weights <- c(w0 = w0, w1 = 0.5 - w3, w2 = 0.5 - w0, w3 = w3)
  errors <- c(w0 = err0, w1 = err3, w2 = err0, w3 = err3)
  standard_errors <- c(w0 = se0, w1 = se3, w2 = se0, w3 = se3)
  if (any(!is.finite(weights)) || any(weights < -1e-7) ||
      any(weights > 1 + 1e-7)) {
    warning("Weights outside [0, 1]; inspect numerical accuracy.", call. = FALSE)
  }
  structure(list(w = weights, abs.error = errors, se = standard_errors,
                 method = method, nsim = nsim), class = "chibar_weights")
}

#' Print chi-bar-square weights
#' @param x An object returned by `chibar_weights()`.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
#' Chi-bar-square weights for a 2-by-2 random-effects covariance test
#'
#' Compute chi-bar-square mixture weights for testing `Gamma = 0` against
#' the full 2-by-2 positive-semidefinite covariance cone. The supplied
#' information must be the **efficient** Fisher information for the ordered
#' coordinates `(gamma11, gamma12, gamma22)`.
#'
#' @param information A symmetric positive-definite 3-by-3 efficient
#'   information matrix, or a `fim_glmm()` result containing `I_metric`.
#' @param method Either `"lorentz"` (solid-angle quadrature) or `"wu_neale"`.
#' @param rel.tol Relative quadrature tolerance.
#' @param subdivisions Positive integer quadrature subdivision limit.
#' @return An object of class `chibar_weights` with components `w`,
#'   `abs.error`, `se`, `method`, and `nsim`.
#' @details The weights correspond to chi-square degrees of freedom 0:3.
#'   The Lorentz method evaluates positive solid angles for the primal and
#'   polar cones. No negative-angle orientation convention is used.
#' @examples
#' chibar_weights(diag(3))
#' @export
print.chibar_weights <- function(x, ...) {
  cat("Chi-bar-square weights (", x$method, ")\n", sep = "")
  print(x$w)
  if (identical(x$method, "counting")) {
    cat("\nMonte Carlo standard errors:\n")
    print(x$se)
    cat("Simulations:", x$nsim, "\n")
  } else {
    cat("\nEstimated absolute integration errors:\n")
    print(x$abs.error)
  }
  invisible(x)
}

.lorentz_transform <- function(information) {
  s <- 1 / sqrt(2)
  T <- rbind(c(s, 0, s), c(s, 0, -s), c(0, sqrt(2), 0))
  .matrix_sqrt_spd(information) %*% solve(T)
}

.lorentz_probability <- function(B, rel.tol, subdivisions) {
  axis <- B[, 1L]
  axis_norm <- sqrt(sum(axis^2))
  a <- axis / axis_norm
  numerator <- abs(det(B)) / axis_norm
  aB <- drop(crossprod(a, B))
  integrand <- function(phi) {
    vapply(phi, function(t) {
      r <- c(1, cos(t), sin(t))
      Br <- drop(B %*% r)
      norm_Br <- sqrt(sum(Br^2))
      aBr <- sum(aB * r)
      numerator / (norm_Br * (norm_Br + aBr))
    }, numeric(1))
  }
  ans <- stats::integrate(integrand, 0, 2 * pi, rel.tol = rel.tol,
                   subdivisions = as.integer(subdivisions),
                   stop.on.error = TRUE)
  list(value = ans$value / (4 * pi), abs.error = ans$abs.error / (4 * pi))
}

.chibar_weights_lorentz <- function(information, rel.tol, subdivisions) {
  B <- .lorentz_transform(information)
  Bpolar <- -solve(t(B))
  primal <- .lorentz_probability(B, rel.tol, subdivisions)
  polar <- .lorentz_probability(Bpolar, rel.tol, subdivisions)
  .assemble_chibar_result(w0 = polar$value, w3 = primal$value,
                         method = "lorentz", err0 = polar$abs.error,
                         err3 = primal$abs.error)
}

.chibar_weights_wu_neale <- function(information, rel.tol, subdivisions) {
  V <- rbind(c(0, 0, 0.5), c(0, -1, 0), c(0.5, 0, 0))
  I_inv_sqrt <- .matrix_sqrt_spd(information, inverse = TRUE)
  H <- I_inv_sqrt %*% V %*% I_inv_sqrt
  ev <- base::eigen((H + t(H)) / 2, symmetric = TRUE)$values
  if (sum(ev > 0) != 1L || sum(ev < 0) != 2L) {
    stop("Unexpected eigenvalue signature.", call. = FALSE)
  }
  lambda <- sort(-ev[ev < 0])
  lambda1 <- lambda[1L]
  lambda2 <- lambda[2L]
  lambda3 <- ev[ev > 0]
  f3 <- function(psi) {
    q <- lambda1 * cos(psi)^2 + lambda2 * sin(psi)^2
    sqrt(q / (lambda3 + q))
  }
  f0 <- function(psi) {
    q <- cos(psi)^2 / lambda1 + sin(psi)^2 / lambda2
    sqrt(q / (1 / lambda3 + q))
  }
  a3 <- stats::integrate(f3, 0, pi / 2, rel.tol = rel.tol,
                  subdivisions = as.integer(subdivisions), stop.on.error = TRUE)
  a0 <- stats::integrate(f0, 0, pi / 2, rel.tol = rel.tol,
                  subdivisions = as.integer(subdivisions), stop.on.error = TRUE)
  .assemble_chibar_result(w0 = 0.5 - a0$value / pi,
                         w3 = 0.5 - a3$value / pi,
                         method = "wu_neale",
                         err0 = a0$abs.error / pi,
                         err3 = a3$abs.error / pi)
}

#' Chi-bar-square weights for a 2-by-2 random-effects covariance test
#'
#' Compute chi-bar-square mixture weights for testing `Gamma = 0` against
#' the full 2-by-2 positive-semidefinite covariance cone. The supplied
#' information must be the **efficient** Fisher information for the ordered
#' coordinates `(gamma11, gamma12, gamma22)`.
#'
#' @param information A symmetric positive-definite 3-by-3 efficient
#'   information matrix, or a `fim_glmm()` result containing `I_metric`.
#' @param method Either `"lorentz"` (solid-angle quadrature) or `"wu_neale"`.
#' @param rel.tol Relative quadrature tolerance.
#' @param subdivisions Positive integer quadrature subdivision limit.
#' @return An object of class `chibar_weights` with components `w`,
#'   `abs.error`, `se`, `method`, and `nsim`.
#' @details The weights correspond to chi-square degrees of freedom 0:3.
#'   The Lorentz method evaluates positive solid angles for the primal and
#'   polar cones. No negative-angle orientation convention is used.
#' @examples
#' chibar_weights(diag(3))
#' @export
chibar_weights <- function(
    information,
    method = c("lorentz", "wu_neale"),
    rel.tol = 1e-9,
    subdivisions = 500L
) {
  
  method <- base::match.arg(method)
  
  # Accept a fim_glmm() result directly
  if (is.list(information) &&
      !is.null(information$I_metric)) {
    information <- information$I_metric
  }
  
  # Validate inputs
  information <- .validate_chibar_information(information)
  
  .validate_chibar_quadrature(
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )
  
  # Select method
  switch(
    method,
    
    lorentz = .chibar_weights_lorentz(
      information = information,
      rel.tol = rel.tol,
      subdivisions = subdivisions
    ),
    
    wu_neale = .chibar_weights_wu_neale(
      information = information,
      rel.tol = rel.tol,
      subdivisions = subdivisions
    )
  )
}
