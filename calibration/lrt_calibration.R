#' Validate chi-bar-square weights
#' @param weights Numeric vector of four mixture weights, ordered df 0:3.
#' @return A normalized numeric vector of weights.
#' @keywords internal
.validate_chibar_weights <- function(weights) {
  if (!is.numeric(weights) || length(weights) != 4L ||
      anyNA(weights) || any(!is.finite(weights)) ||
      any(weights < -1e-8) || abs(sum(weights) - 1) > 1e-6) {
    stop("weights must be four finite, nonnegative (up to tolerance) values summing to one.", call. = FALSE)
  }
  weights <- pmax(weights, 0)
  weights / sum(weights)
}

#' Chi-bar-square upper-tail probability
#'
#' Evaluates P(X >= t) for t = 0 and P(X > t) for t > 0,
#' where X has a chi-bar-square distribution with df 0, 1, 2, 3.
#' This convention returns one at t = 0, as in the original LRT study.
#'
#' @param statistic Numeric vector of nonnegative likelihood-ratio statistics.
#' @param weights Four mixture weights in degrees-of-freedom order 0:3.
#' @return Numeric vector of upper-tail probabilities.
#' @export
chibar_survival <- function(statistic, weights) {
  w <- .validate_chibar_weights(weights)
  if (!is.numeric(statistic) || anyNA(statistic) ||
      any(!is.finite(statistic)) || any(statistic < 0)) {
    stop("statistic must contain finite nonnegative values.", call. = FALSE)
  }
  vapply(statistic, function(t) {
    if (t == 0) return(1)
    sum(w[2:4] * stats::pchisq(t, df = 1:3, lower.tail = FALSE))
  }, numeric(1))
}

#' Chi-bar-square critical value
#'
#' Returns the nonnegative cutoff c such that rejection for T > c
#' has probability no larger than alpha under the reference mixture.
#' If alpha >= 1 - w0, the cutoff is zero and the test is conservative.
#'
#' @param alpha Scalar nominal level in (0,1).
#' @param weights Four chi-bar-square weights, df 0:3.
#' @return Nonnegative critical value.
#' @export
chibar_critical_value <- function(alpha, weights) {
  w <- .validate_chibar_weights(weights)
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0,1).", call. = FALSE)
  }
  if (alpha >= 1 - w[1L]) return(0)
  f <- function(t) sum(w[2:4] * stats::pchisq(t, df = 1:3, lower.tail = FALSE)) - alpha
  upper <- 1
  while (f(upper) > 0) {
    upper <- upper * 2
    if (!is.finite(upper)) {
      stop("Unable to bracket chi-bar-square critical value.", call. = FALSE)
    }
  }
  stats::uniroot(f, c(0, upper), tol = 1e-10)$root
}

#' Calibrate a global covariance-component LRT
#'
#' @param statistic Nonnegative LRT statistic(s).
#' @param weights Chi-bar-square weights ordered df 0:3.
#' @param alpha Nominal significance level(s) in (0,1).
#' @return Data frame with statistic, alpha, p_value, critical_value and reject.
#' @export
calibrate_global_lrt <- function(statistic, weights, alpha = 0.05) {
  if (!is.numeric(alpha) || !length(alpha) || anyNA(alpha) ||
      any(!is.finite(alpha)) || any(alpha <= 0 | alpha >= 1)) {
    stop("alpha must contain values in (0,1).", call. = FALSE)
  }
  if (!is.numeric(statistic) || !length(statistic) || anyNA(statistic) ||
      any(!is.finite(statistic)) || any(statistic < 0)) {
    stop("statistic must contain finite nonnegative values.", call. = FALSE)
  }
  p <- chibar_survival(statistic, weights)
  cutoffs <- vapply(alpha, chibar_critical_value, numeric(1), weights = weights)
  out <- expand.grid(statistic = statistic, alpha = alpha, KEEP.OUT.ATTRS = FALSE)
  out$p_value <- rep(p, times = length(alpha))
  out$critical_value <- rep(cutoffs, each = length(statistic))
  # A nonrandomized LRT with an atom at zero cannot reject at T = 0.
  out$reject <- out$statistic > out$critical_value
  out
}
