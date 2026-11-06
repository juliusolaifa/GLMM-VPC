# Load the consolidated fitter, revised likelihood adapter, and vertex adapter.
testthat::test_that('pointwise closed-null wrapper validates arguments', {
  testthat::expect_error(fit_vpc_pointwise_null(list(), 0), 'glmmVPC_fit')
  fake <- structure(list(engine='glmmTMB', raw_fit=list(), control=list()), class='glmmVPC_fit')
  testthat::expect_error(fit_vpc_pointwise_null(fake, Inf), 'finite')
  testthat::expect_error(fit_vpc_pointwise_null(fake, 0, likelihood_tol=-1), 'likelihood_tol')
})

testthat::test_that('pointwise closed-null compares ray and vertex likelihoods', {
  testthat::skip_if_not_installed('glmmTMB')
  testthat::skip_if_not_installed('lme4')
  set.seed(20261105)
  n <- 30L; m <- 8L
  id <- factor(rep(seq_len(n), each=m))
  x <- rep(rep(c(0,1), each=m/2), n)
  u <- rnorm(n, sd=.5)
  y <- rpois(length(x), exp(.2+.3*x+u[as.integer(id)]))
  dat <- data.frame(y=y,x=x,id=id)
  fit <- fit_glmm(y ~ x + (1+x|id), data=dat, family='poisson')
  ray <- tryCatch(fit_vpc_pointwise_ray(fit, 0, FALSE), error=function(e)e)
  vertex <- tryCatch(fit_vpc_omnibus_null(fit, FALSE), error=function(e)e)
  testthat::skip_if(inherits(ray,'error') || inherits(vertex,'error'),
                    'Restricted optimizer unavailable on this platform')
  null <- fit_vpc_pointwise_null(fit, 0, require_acceptable=FALSE)
  testthat::expect_equal(null$logLik, max(ray$logLik,vertex$logLik), tolerance=1e-6)
  testthat::expect_true(null$diagnostics$boundary_complete)
  testthat::expect_false(null$diagnostics$global_optimality_verified)
  testthat::expect_true(null$diagnostics$selected_stratum %in% c('ray','vertex'))
  testthat::expect_equal(drop(crossprod(c(1,0), null$Gamma %*% c(1,0))),
                         0, tolerance=1e-8)
})
