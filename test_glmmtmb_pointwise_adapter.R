# Source fit_glmm_consolidated.R, constrained_fit_revised.R,
# glmmtmb_likelihood_adapter_revised.R, and
# glmmtmb_pointwise_adapter_revised.R before running standalone.

testthat::test_that('pointwise adapter rejects invalid inputs', {
  testthat::expect_error(fit_vpc_pointwise_null(list(), 0), 'glmmVPC_fit')
  fake <- structure(list(engine='glmmTMB'), class='glmmVPC_fit')
  testthat::expect_error(fit_vpc_pointwise_null(fake, Inf), 'finite')
  testthat::expect_error(fit_vpc_pointwise_null(fake, 0, NA),
                         'require_acceptable')
})

testthat::test_that('pointwise ray has exact zero variance at x0', {
  for (x0 in c(-2, 0, .5, 1, 3)) {
    direction <- c(-x0, 1)
    Gamma <- .7 * tcrossprod(direction)
    testthat::expect_equal(drop(crossprod(c(1,x0), Gamma %*% c(1,x0))),
                           0, tolerance=1e-12)
    testthat::expect_equal(det(Gamma), 0, tolerance=1e-12)
  }
})

testthat::test_that('pointwise adapter preserves rows and flags vertex omission', {
  testthat::skip_if_not_installed('glmmTMB')
  testthat::skip_if_not_installed('lme4')
  set.seed(20261105)
  n <- 50L; m <- 10L
  id <- factor(rep(seq_len(n), each=m))
  x <- rep(rep(c(0,1), each=m/2), n)
  u <- rnorm(n, sd=.5)
  y <- rpois(length(x), exp(.2+.3*x+u[as.integer(id)]))
  dat <- data.frame(y=y, x=x, id=id)
  fit <- fit_glmm(y ~ x + (1+x|id), data=dat, family='poisson')
  for (x0 in c(0, 1, .5)) {
    null <- fit_vpc_pointwise_null(fit, x0, require_acceptable=FALSE)
    z <- c(1,x0)
    testthat::expect_equal(drop(crossprod(z, null$Gamma %*% z)), 0,
                           tolerance=1e-8)
    testthat::expect_true(null$diagnostics$observation_rows_identical)
    testthat::expect_true(null$diagnostics$native_objective_consistency$consistent)
    testthat::expect_false(null$diagnostics$boundary_complete)
    testthat::expect_false(null$diagnostics$global_optimality_verified)
    testthat::expect_equal(null$objective, -null$logLik, tolerance=1e-4)
    testthat::expect_true(is.finite(null$logLik))
  }
})
