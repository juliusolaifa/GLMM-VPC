#' Variance Partition Coefficient for a Random-Intercept-and-Slope GLMM
#'
#' Computes the observation-scale variance partition coefficient at each
#' supplied covariate value for a log-link count GLMM or identity-link
#' Gaussian mixed model.
#'
#' @param beta Numeric vector of length two: fixed intercept and slope.
#' @param Gamma Symmetric positive-semidefinite 2-by-2 random-effects
#'   covariance matrix, ordered intercept then slope.
#' @param x Numeric vector of covariate values.
#' @param family One of `"poisson"`, `"nbinom"`, `"cpg"`, or `"gaussian"`.
#' @param family_params Named numeric vector: empty for Poisson,
#'   `c(theta = ...)` for negative binomial,
#'   `c(phi = ..., power = ...)` for compound Poisson-Gamma, or
#'   `c(sigma2 = ...)` for Gaussian.
#'
#' @return Numeric vector of VPC values, in the same order as `x`.
#' @export
#' @examples
#' glmm_vpc(c(0.5, 0.2), diag(c(0.2, 0.1)), c(0, 1), "poisson")
glmm_vpc <- function(
    beta,
    Gamma,
    x,
    family = c("poisson", "nbinom", "cpg", "gaussian"),
    family_params = numeric(0)
) {
  family <- match.arg(family)

  if (!is.numeric(beta) || length(beta) != 2L ||
      any(!is.finite(beta))) {
    stop("`beta` must be a finite numeric vector of length two.", call. = FALSE)
  }
  if (!is.matrix(Gamma) || !is.numeric(Gamma) ||
      !identical(dim(Gamma), c(2L, 2L)) || any(!is.finite(Gamma))) {
    stop("`Gamma` must be a finite numeric 2-by-2 matrix.", call. = FALSE)
  }
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop("`x` must be a finite numeric vector.", call. = FALSE)
  }
  if (!is.numeric(family_params) || any(!is.finite(family_params))) {
    stop("`family_params` must be a finite numeric vector.", call. = FALSE)
  }
  if (max(abs(Gamma - t(Gamma))) > 1e-8) {
    stop("`Gamma` must be symmetric.", call. = FALSE)
  }
  if (min(eigen(Gamma, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
    stop("`Gamma` must be positive semidefinite.", call. = FALSE)
  }

  eta <- beta[1L] + beta[2L] * x
  v <- Gamma[1L, 1L] + 2 * x * Gamma[1L, 2L] +
    x^2 * Gamma[2L, 2L]
  v <- pmax(v, 0)

  log_A <- rep(-Inf, length(v))
  positive <- v > 0
  log_A[positive] <- log(expm1(v[positive]))

  if (family == "gaussian") {
    if (!identical(names(family_params), "sigma2") ||
        family_params[1L] <= 0) {
      stop("Gaussian requires `family_params = c(sigma2 = ...)` with sigma2 > 0.",
           call. = FALSE)
    }
    log_between <- log(v)
    log_within <- rep(log(unname(family_params["sigma2"])), length(x))
  } else {
    log_between <- 2 * eta + v + log_A

    if (family == "poisson") {
      if (length(family_params) != 0L) {
        stop("Poisson requires no family parameters.", call. = FALSE)
      }
      log_within <- eta + v / 2
    } else if (family == "nbinom") {
      if (!identical(names(family_params), "theta") ||
          family_params[1L] <= 0) {
        stop("Negative binomial requires `family_params = c(theta = ...)` with theta > 0.",
             call. = FALSE)
      }
      theta <- unname(family_params["theta"])
      a <- eta + v / 2
      b <- 2 * eta + 2 * v - log(theta)
      log_within <- pmax(a, b) + log1p(exp(-abs(a - b)))
    } else if (family == "cpg") {
      if (length(family_params) != 2L ||
          !setequal(names(family_params), c("phi", "power"))) {
        stop("CPG requires exactly `phi` and `power`.", call. = FALSE)
      }
      phi <- unname(family_params["phi"])
      power <- unname(family_params["power"])
      if (phi <= 0 || power <= 1 || power >= 2) {
        stop("CPG requires phi > 0 and 1 < power < 2.", call. = FALSE)
      }
      log_within <- log(phi) + power * eta + 0.5 * power^2 * v
    }
  }

  stats::plogis(log_between - log_within)
}
