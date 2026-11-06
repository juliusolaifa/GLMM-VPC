test_that("chi-bar survival respects the atom at zero", {
  w <- c(.2, .3, .3, .2)
  expect_equal(chibar_survival(0, w), 1)
  expect_equal(chibar_survival(1e-8, w),
               sum(w[2:4] * pchisq(1e-8, df=1:3, lower.tail=FALSE)))
})

test_that("critical values attain the nominal positive-tail level", {
  w <- c(.2, .3, .3, .2)
  for (alpha in c(.01, .05, .1, .5)) {
    cutoff <- chibar_critical_value(alpha, w)
    expect_equal(chibar_survival(cutoff, w), alpha, tolerance=1e-8)
  }
  expect_equal(chibar_critical_value(.85, w), 0)
})

test_that("calibration is conservative at the atom and rejects strictly", {
  w <- c(.2, .3, .3, .2)
  out <- calibrate_global_lrt(c(0, 1, 10), w, alpha=c(.05,.85))
  expect_equal(nrow(out), 6L)
  expect_false(any(out$reject[out$statistic == 0]))
  expect_error(calibrate_global_lrt(numeric(), w), "statistic")
  expect_error(chibar_survival(-1, w), "statistic")
  expect_error(chibar_survival(Inf, w), "statistic")
})
