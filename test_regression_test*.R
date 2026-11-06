# Run after sourcing fit_glmm_consolidated.R, constrained_fit_revised.R,
# and glmmtmb_likelihood_adapter_revised.R (or in the package test environment).

testthat::test_that("native likelihood rejects malformed fit and parameters", {
  testthat::expect_error(vpc_glmmtmb_likelihood(list()), "glmmVPC_fit")
  fake <- structure(list(engine="glmmTMB",raw_fit=list(),control=list(ll_tol=1e-6)),
                    class="glmmVPC_fit")
  testthat::expect_error(vpc_glmmtmb_likelihood(fake), "does not expose")
})

testthat::test_that("row extraction preserves analysis row order", {
  a <- list(frame=data.frame(y=1:3, row.names=c("1","3","4")))
  testthat::expect_identical(.vpc_glmmtmb_rows(a), c("1","3","4"))
  testthat::expect_null(.vpc_glmmtmb_rows(list()))
})

testthat::test_that("omnibus null uses exact zero covariance", {
  testthat::skip_if_not_installed("glmmTMB")
  testthat::skip_if_not_installed("lme4")
  set.seed(20261105)
  id <- factor(rep(seq_len(40), each=10))
  x <- rep(rep(c(0,1), each=5), 40)
  u <- rnorm(40, sd=.35)
  y <- rpois(length(x), exp(.3 + .25*x + u[as.integer(id)]))
  dat <- data.frame(y=y,x=x,id=id)
  fit <- fit_glmm(y ~ x + (1+x|id), dat, family="poisson",
                  retry_mode="strict")
  lik <- vpc_glmmtmb_likelihood(fit)
  testthat::expect_equal(lik$nll(lik$start), -fit$logLik,
                         tolerance=1e-4)
  testthat::expect_length(lik$gradient(lik$start), length(lik$start))
  null <- fit_vpc_omnibus_null(fit, require_acceptable=FALSE)
  testthat::expect_true(all(null$Gamma == 0))
  testthat::expect_true(null$diagnostics$observation_rows_identical)
  testthat::expect_true(null$diagnostics$boundary_complete)
  testthat::expect_false(null$diagnostics$global_optimality_verified)
  testthat::expect_equal(null$objective, -null$logLik,
                         tolerance=1e-4)
  testthat::expect_gte(fit$logLik - null$logLik, -1e-6)
})
