# Run after sourcing constrained_fit_revised.R, or within a package test suite.
testthat::test_that("exact map and independent numeric validation", {
  p <- vpc_constraint(
    nll=function(z) (z[1]-2)^2+(z[2]-4)^2,
    map=function(u) c(u,2*u),start=0,
    validate=function(z) abs(z[2]-2*z[1]),label="linear equality")
  f <- fit_glmm_constrained(p,require_validation=TRUE,require_acceptable=TRUE)
  testthat::expect_equal(f$par,c(2,4),tolerance=1e-3)
  testthat::expect_true(f$diagnostics$feasible)
  testthat::expect_true(f$diagnostics$stationary)
  testthat::expect_equal(f$diagnostics$constraint_violation,0,tolerance=1e-8)
  testthat::expect_false(f$diagnostics$boundary_complete)
})

testthat::test_that("nonlinear map uses chain rule for gradient diagnostics", {
  p <- vpc_constraint(function(z) (z-4)^2, function(u) u^2, start=1,
                      validate=function(z) 0)
  f <- fit_glmm_constrained(p,grad=function(z) 2*(z-4),
                            require_validation=TRUE,require_acceptable=TRUE)
  testthat::expect_equal(f$free_par,2,tolerance=1e-3)
  testthat::expect_true(f$diagnostics$stationary)
})

testthat::test_that("missing validation does not certify numerical acceptability", {
  p <- vpc_constraint(function(z) (z-1)^2,identity,start=0)
  f <- fit_glmm_constrained(p)
  testthat::expect_true(is.na(f$diagnostics$feasible))
  testthat::expect_false(f$diagnostics$numerically_acceptable)
  testthat::expect_error(fit_glmm_constrained(p,require_validation=TRUE),
                         "validate")
})

testthat::test_that("feasibility tolerance applies to numeric validator", {
  p <- vpc_constraint(function(z) z^2,identity,start=1,
                      validate=function(z) 0.01)
  testthat::expect_error(fit_glmm_constrained(p,feasibility_tol=1e-4),
                         "No finite feasible")
  f <- fit_glmm_constrained(p,feasibility_tol=0.02)
  testthat::expect_true(f$diagnostics$feasible)
})

testthat::test_that("LRT protects against negative values", {
  testthat::expect_equal(vpc_lrt_statistic(list(logLik=-10),list(logLik=-11)),2)
  testthat::expect_warning(
    testthat::expect_equal(vpc_lrt_statistic(list(logLik=-10),
                                           list(logLik=-9.9999999)),0),"Small negative")
  testthat::expect_error(vpc_lrt_statistic(list(logLik=-10),list(logLik=-9)),
                         "Restricted log-likelihood")
})

testthat::test_that("zero-dimensional optimization is supported", {
  p <- vpc_constraint(function(z) sum(z^2),function(u) c(0,0),
                      start=numeric(),validate=function(z) all(z==0))
  f <- fit_glmm_constrained(p,require_validation=TRUE,require_acceptable=TRUE)
  testthat::expect_equal(f$objective,0)
  testthat::expect_equal(f$diagnostics$max_abs_reduced_gradient,0)
})
