simulate_s_dist <- function(model_fn, ..., n_blocks = 20000, kmax = 80) {
  p <- model_fn(..., k = 0:kmax)
  p <- p / sum(p)
  counts <- as.vector(stats::rmultinom(1, n_blocks, p))
  data.frame(k = 0:kmax, count = counts)
}

test_that("fit_im recovers known parameters from simulated data", {
  set.seed(42)
  s_dist <- simulate_s_dist(wh_im, M = 2.0, tau = 1.2, theta = 3.0)
  fit <- fit_im(s_dist)
  expect_equal(unname(fit$par["M"]), 2.0, tolerance = 0.2)
  expect_equal(unname(fit$par["tau"]), 1.2, tolerance = 0.1)
  expect_equal(unname(fit$par["theta"]), 3.0, tolerance = 0.1)
  expect_equal(fit$convergence, 0)
})

test_that("fit_si recovers known parameters from simulated data", {
  set.seed(43)
  s_dist <- simulate_s_dist(wh_si, tau = 0.6, theta = 1.5)
  fit <- fit_si(s_dist)
  expect_equal(unname(fit$par["tau"]), 0.6, tolerance = 0.1)
  expect_equal(unname(fit$par["theta"]), 1.5, tolerance = 0.1)
})

test_that("user-supplied bounds override defaults for a single parameter", {
  s_dist <- data.frame(k = 0:5, count = c(100, 45, 20, 9, 4, 2))
  fit_default <- fit_im(s_dist)
  fit_narrow <- fit_im(s_dist, bounds = list(M = c(0.001, 0.01)))
  expect_true(fit_narrow$par[["M"]] <= 0.01)
  expect_true(fit_narrow$bounds$M[2] == 0.01)
  expect_true(fit_default$bounds$M[2] == default_bounds$IM$M[2])
})

test_that("unknown bound names raise an informative error", {
  s_dist <- data.frame(k = 0:5, count = c(100, 45, 20, 9, 4, 2))
  expect_error(fit_im(s_dist, bounds = list(not_a_param = c(0, 1))), "Unknown bound")
})

test_that("S3 methods coef/logLik/print/summary work as expected", {
  s_dist <- data.frame(k = 0:5, count = c(100, 45, 20, 9, 4, 2))
  fit <- fit_im(s_dist)

  expect_equal(coef(fit), fit$par)
  expect_equal(logLik(fit), fit$loglik)
  expect_false(inherits(logLik(fit), "logLik"))  # deliberately NOT AIC/BIC-compatible

  expect_output(print(fit), "Isolation-with-Migration model")
  expect_output(print(summary(fit)), "Call:")
  expect_output(print(summary(fit)), "Parameter estimates and bounds")
})

test_that("plot.demogfit_fit runs without error", {
  s_dist <- data.frame(k = 0:5, count = c(100, 45, 20, 9, 4, 2))
  fit <- fit_im(s_dist)
  pdf(NULL)  # avoid opening a graphics device / writing a file
  on.exit(dev.off())
  expect_no_error(plot(fit))
})

test_that("this package does not provide AIC()/BIC() methods", {
  expect_false(exists("AIC.demogfit_fit"))
  expect_false(exists("BIC.demogfit_fit"))
})

test_that("new default bounds: tau-type parameters floored at 0.5, M ceiled at 5, f untouched", {
  expect_equal(default_bounds$SI$tau[1], 0.5)
  expect_equal(default_bounds$IM$tau[1], 0.5)
  expect_equal(default_bounds$IM$M[2], 5)
  expect_equal(default_bounds$IIM$tau1[1], 0.5)
  expect_equal(default_bounds$IIM$tau0[1], 0.5)
  expect_equal(default_bounds$IIM$M[2], 5)
  expect_equal(default_bounds$SC$tau1[1], 0.5)
  expect_equal(default_bounds$SC$tau0[1], 0.5)
  expect_equal(default_bounds$SC$f, c(1e-3, 1))       # untouched - not on M's scale
})

test_that("parameters shared by SI/IM/IIM have IDENTICAL bounds, so SI-vs-IM and IM-vs-IIM LRTs are not biased by one model having a narrower search range than the other", {
  # tau: shared by SI and IM
  expect_equal(default_bounds$SI$tau, default_bounds$IM$tau)
  # theta: shared by all three
  expect_equal(default_bounds$SI$theta, default_bounds$IM$theta)
  expect_equal(default_bounds$IM$theta, default_bounds$IIM$theta)
  # M: shared by IM and IIM
  expect_equal(default_bounds$IM$M, default_bounds$IIM$M)
  # SC is NOT nested with SI/IM/IIM (see fit_demography), so it is
  # deliberately exempt from this consistency requirement - confirm it's
  # untouched by this round of changes rather than assert it matches
  expect_equal(default_bounds$SC$f, c(1e-3, 1))
  expect_equal(default_bounds$SC$theta, c(1e-3, 8))
})

test_that("boundary_hits is empty for a normal, well-converged interior fit", {
  set.seed(20)
  s_dist <- simulate_s_dist(wh_im, M = 1.5, tau = 3, theta = 2, n_blocks = 8000)
  fit <- fit_im(s_dist)
  expect_length(fit$boundary_hits, 0)
  expect_false(any(grepl("WARNING", capture.output(print(fit)))))
})

test_that("boundary_hits correctly detects a parameter forced to its bound", {
  set.seed(21)
  s_dist <- simulate_s_dist(wh_im, M = 1.5, tau = 3, theta = 2, n_blocks = 8000)
  fit <- fit_im(s_dist, bounds = list(tau = c(2.9, 3.0)))
  expect_true("tau" %in% names(fit$boundary_hits))
  expect_true(fit$boundary_hits[["tau"]] %in% c("lower", "upper"))
  out <- capture.output(print(fit))
  expect_true(any(grepl("WARNING", out)))
  out_summary <- capture.output(print(summary(fit)))
  expect_true(any(grepl("AT BOUND", out_summary)))
})

test_that("boundary_hits respects user-overridden bounds, not just defaults", {
  set.seed(22)
  s_dist <- simulate_s_dist(wh_im, M = 1.5, tau = 3, theta = 2, n_blocks = 8000)
  # tau=3 is comfortably interior to the DEFAULT bound (0.5, 15) - should
  # only be flagged once the bound is narrowed to make it boundary-adjacent
  fit_default <- fit_im(s_dist)
  expect_false("tau" %in% names(fit_default$boundary_hits))
  fit_narrow <- fit_im(s_dist, bounds = list(tau = c(2.95, 3.05)))
  expect_true("tau" %in% names(fit_narrow$boundary_hits))
})

test_that("fit_demography()'s print/summary surface boundary hits found in ANY fitted model", {
  set.seed(23)
  s_dist <- simulate_s_dist(wh_im, M = 1.5, tau = 3, theta = 2, n_blocks = 8000)
  res <- fit_demography(s_dist, models = c("SI", "IM"),
                         bounds = list(IM = list(tau = c(2.9, 3.0))))
  out <- capture.output(print(res))
  expect_true(any(grepl("WARNING", out)))
  expect_true(any(grepl("IM:.*tau", out)))
  out_s <- capture.output(print(summary(res)))
  expect_true(any(grepl("WARNING", out_s)))
})
