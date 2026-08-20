sim_pair <- function(M, tau, theta, n = 8000, kmax = 60) {
  p <- wh_im(M = M, tau = tau, theta = theta, k = 0:kmax)
  data.frame(k = 0:kmax, count = as.vector(stats::rmultinom(1, n, p / sum(p))))
}

test_that("scale_parameters selects SI when there is no evidence for gene flow", {
  set.seed(1)
  s_dist <- sim_pair(M = 0, tau = 0.7, theta = 2, n = 6000)
  result <- fit_demography(s_dist)
  scaled <- scale_parameters(result, mu = 2.8e-9, block_length = 200)
  expect_equal(scaled$model, "SI")
  expect_equal(scaled$migration_rate_per_generation, 0)
  expect_true(scaled$Ne > 0)
  expect_true(scaled$divergence_time_generations > 0)
})

test_that("scale_parameters selects IM when there is gene flow but no IIM support", {
  set.seed(2)
  s_dist <- sim_pair(M = 1.5, tau = 1, theta = 2, n = 6000)
  result <- fit_demography(s_dist, models = c("SI", "IM"))  # IIM not fit at all
  scaled <- scale_parameters(result, mu = 2.8e-9, block_length = 200)
  expect_equal(scaled$model, "IM")
  expect_true(scaled$migration_rate_per_generation > 0)
  expect_true(is.na(scaled$migration_duration_generations))
})

test_that("scale_parameters never selects SC", {
  set.seed(3)
  s_dist <- sim_pair(M = 2, tau = 1, theta = 2, n = 10000)
  result <- fit_demography(s_dist)  # all four models fit
  scaled <- scale_parameters(result, mu = 2.8e-9, block_length = 200)
  expect_true(scaled$model %in% c("SI", "IM", "IIM"))
})

test_that("generation_time_years is optional and converts correctly", {
  set.seed(4)
  s_dist <- sim_pair(M = 1.5, tau = 1, theta = 2, n = 6000)
  result <- fit_demography(s_dist, models = c("SI", "IM"))
  scaled_no_gt <- scale_parameters(result, mu = 2.8e-9, block_length = 200)
  scaled_gt <- scale_parameters(result, mu = 2.8e-9, block_length = 200, generation_time_years = 0.1)
  expect_null(scaled_no_gt$divergence_time_years)
  expect_false(is.null(scaled_gt$divergence_time_years))
  expect_equal(scaled_gt$divergence_time_years, scaled_gt$divergence_time_generations * 0.1)
})

test_that("ploidy_scalar affects Ne but not divergence time or migration rate", {
  set.seed(5)
  s_dist <- sim_pair(M = 1.5, tau = 1, theta = 2, n = 6000)
  result <- fit_demography(s_dist, models = c("SI", "IM"))
  s4 <- scale_parameters(result, mu = 2.8e-9, block_length = 200, ploidy_scalar = 4)
  s2 <- scale_parameters(result, mu = 2.8e-9, block_length = 200, ploidy_scalar = 2)
  expect_equal(s2$Ne, 2 * s4$Ne)  # halving the scalar doubles Ne
  expect_equal(s4$divergence_time_generations, s2$divergence_time_generations)
  expect_equal(s4$migration_rate_per_generation, s2$migration_rate_per_generation)
})

test_that("scale_parameters requires mu and block_length explicitly", {
  set.seed(6)
  s_dist <- sim_pair(M = 1, tau = 1, theta = 2, n = 4000)
  result <- fit_demography(s_dist, models = c("SI", "IM"))
  expect_error(scale_parameters(result, block_length = 200), "must both be supplied")
  expect_error(scale_parameters(result, mu = 2.8e-9), "must both be supplied")
})

test_that("scale_parameters errors on non-demogfit_result input", {
  expect_error(scale_parameters(list(), mu = 1e-9, block_length = 100), "demogfit_result")
})

test_that("print.demogfit_scaled runs without error", {
  set.seed(7)
  s_dist <- sim_pair(M = 1.5, tau = 1, theta = 2, n = 6000)
  result <- fit_demography(s_dist, models = c("SI", "IM"))
  scaled <- scale_parameters(result, mu = 2.8e-9, block_length = 200, generation_time_years = 0.1)
  expect_output(print(scaled), "Scaled parameter estimates")
  expect_output(print(scaled), "Effective population size")
})
