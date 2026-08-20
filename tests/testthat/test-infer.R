sim_pair <- function(M, tau, theta, n = 8000, kmax = 60) {
  p <- wh_im(M = M, tau = tau, theta = theta, k = 0:kmax)
  data.frame(k = 0:kmax, count = as.vector(stats::rmultinom(1, n, p / sum(p))))
}

test_that("fit_demography detects gene flow when data simulated under IM", {
  set.seed(1)
  s_dist <- sim_pair(M = 2.0, tau = 1.0, theta = 2.5)
  res <- fit_demography(s_dist)
  expect_true(res$gene_flow)
  expect_true(res$gene_flow_test$p_value < 0.05)
  expect_match(res$gene_flow_test$interpretation, "Yes")
})

test_that("fit_demography detects no gene flow when data simulated under SI", {
  set.seed(2)
  s_dist <- sim_pair(M = 0, tau = 0.8, theta = 2.5)
  res <- fit_demography(s_dist)
  expect_false(res$gene_flow)
  expect_match(res$gene_flow_test$interpretation, "No")
})

test_that("fit_demography works with only SI and IM requested", {
  set.seed(3)
  s_dist <- sim_pair(M = 4, tau = 1.0, theta = 2, n = 15000)
  res <- fit_demography(s_dist, models = c("SI", "IM"))
  expect_true(res$gene_flow)
  expect_null(res$iim_test)
  expect_null(res$fits$SC)
})

test_that("fit_demography never statistically compares SC to IIM/IM/SI", {
  set.seed(5)
  s_dist <- sim_pair(M = 2, tau = 1, theta = 2, n = 10000)
  res <- fit_demography(s_dist)
  expect_true(is.list(res$fits$SC))
  # SC must not appear inside either LRT test's model table
  expect_false("SC" %in% res$gene_flow_test$table$model)
  expect_false("SC" %in% res$iim_test$table$model)
  # there should be no p-value/LRT field anywhere named after SC
  expect_null(res$sc_test)
})

test_that("fit_demography errors if SI or IM is missing from models", {
  s_dist <- sim_pair(M = 1, tau = 1, theta = 2)
  expect_error(fit_demography(s_dist, models = c("IIM", "SC")), "at least")
})

test_that("fit_demography_all batches multiple pairs into one table with full results retained", {
  set.seed(4)
  pair_data <- list(
    a = sim_pair(M = 2, tau = 1, theta = 2),
    b = sim_pair(M = 0, tau = 0.7, theta = 2)
  )
  res <- fit_demography_all(pair_data, models = c("SI", "IM"))
  expect_equal(nrow(res), 2)
  expect_setequal(res$pair, c("a", "b"))
  full <- attr(res, "results")
  expect_setequal(names(full), c("a", "b"))
  expect_s3_class(full$a, "demogfit_result")
})

test_that("print and summary methods run without error", {
  set.seed(6)
  s_dist <- sim_pair(M = 2, tau = 1, theta = 2, n = 6000)
  res <- fit_demography(s_dist)
  expect_output(print(res), "Gene flow")
  expect_output(print(summary(res)), "Demographic inference summary")
})

test_that("large block-count caveat is surfaced", {
  set.seed(7)
  s_dist <- sim_pair(M = 2, tau = 1, theta = 2, n = 5000)
  res <- fit_demography(s_dist, no_blocks = 100000)
  expect_output(print(res), "Note: this pair has 100000 blocks")
})

sim_sc_pair <- function(f, tau1, tau0, theta, n = 8000, kmax = 80) {
  p <- wh_sc(f = f, tau1 = tau1, tau0 = tau0, theta = theta, k = 0:kmax)
  data.frame(k = 0:kmax, count = as.vector(stats::rmultinom(1, n, p / sum(p))))
}

test_that("comparison_method defaults to 'lrt' and leaves existing fields untouched", {
  set.seed(10)
  s_dist <- sim_pair(M = 1.5, tau = 1, theta = 2)
  res <- fit_demography(s_dist)
  expect_null(res$model_comparison)
  expect_equal(res$comparison_method, "lrt")

  set.seed(10)
  res_matched <- fit_demography(s_dist, comparison_method = "lrt")
  set.seed(10)
  s_dist2 <- sim_pair(M = 1.5, tau = 1, theta = 2)
  set.seed(10)
  res2 <- fit_demography(s_dist2)
  expect_identical(res_matched$gene_flow_test, res2$gene_flow_test)
})

test_that("comparison_method='aic' ranks every fitted model including SC", {
  set.seed(11)
  s_dist <- sim_pair(M = 1.5, tau = 1, theta = 2)
  res <- fit_demography(s_dist, comparison_method = "aic")
  expect_false(is.null(res$model_comparison))
  expect_setequal(res$model_comparison$table$model, c("SI", "IM", "IIM", "SC"))
  expect_true(res$model_comparison$best_model %in% c("SI", "IM", "IIM", "SC"))
  # AIC = -2*logLik + 2*k, spot-checked for one row
  row <- res$model_comparison$table[res$model_comparison$table$model == "SI", ]
  expect_equal(row$AIC, -2 * row$logLik + 2 * row$k)
})

test_that("AIC correctly identifies SC as best model on data simulated under SC", {
  set.seed(12)
  s_dist <- sim_sc_pair(f = 0.3, tau1 = 0.5, tau0 = 1.2, theta = 2, n = 10000)
  res <- fit_demography(s_dist, comparison_method = "aic")
  expect_equal(res$model_comparison$best_model, "SC")
})

test_that("comparison_method='loglik' flags a caveat only when parameter counts differ", {
  set.seed(13)
  s_dist <- sim_pair(M = 1.5, tau = 1, theta = 2)
  res_all <- fit_demography(s_dist, comparison_method = "loglik")
  expect_false(is.null(res_all$model_comparison$caveat))

  res_iim_sc <- fit_demography(s_dist, models = c("SI", "IM", "IIM", "SC"))
  ranked <- demogfit:::.rank_models(res_iim_sc$fits[c("IIM", "SC")], "loglik")
  expect_null(ranked$caveat)
})

test_that("AIC and raw log-likelihood rankings agree when parameter counts are equal", {
  set.seed(14)
  s_dist <- sim_pair(M = 1.5, tau = 1, theta = 2)
  res <- fit_demography(s_dist, comparison_method = "aic")
  fits_iim_sc <- res$fits[c("IIM", "SC")]
  rank_aic <- demogfit:::.rank_models(fits_iim_sc, "aic")
  rank_loglik <- demogfit:::.rank_models(fits_iim_sc, "loglik")
  expect_identical(rank_aic$best_model, rank_loglik$best_model)
})

test_that("fit_demography_all surfaces best_model when comparison_method is set", {
  set.seed(15)
  pair_data <- list(a = sim_pair(M = 1.5, tau = 1, theta = 2))
  res <- fit_demography_all(pair_data, comparison_method = "aic")
  expect_true("best_model" %in% names(res))
})

test_that("print methods include the model comparison section when present", {
  set.seed(16)
  s_dist <- sim_pair(M = 1.5, tau = 1, theta = 2, n = 6000)
  res <- fit_demography(s_dist, comparison_method = "aic")
  expect_output(print(res), "Model comparison")
  expect_output(print(summary(res)), "Model comparison")
})
