k <- 0:80

test_that("each model's likelihood sums to ~1 over k", {
  expect_equal(sum(wh_si(tau = 1.0, theta = 2.0, k = k)), 1, tolerance = 1e-4)
  expect_equal(sum(wh_im(M = 1.5, tau = 0.8, theta = 2.0, k = k)), 1, tolerance = 1e-4)
  expect_equal(sum(wh_iim(M = 1.2, tau1 = 0.3, tau0_total = 2.0, theta = 2.0, k = k)), 1, tolerance = 1e-4)
  expect_equal(sum(wh_sc(f = 0.3, tau1 = 0.2, tau0 = 1.5, theta = 2.0, k = k)), 1, tolerance = 1e-4)
})

test_that("SI is exactly the M = 0 special case of IM", {
  p_si <- wh_si(tau = 1.0, theta = 2.0, k = k)
  p_im0 <- wh_im(M = 0, tau = 1.0, theta = 2.0, k = k)
  expect_equal(p_si, p_im0)
})

test_that("negative M returns all zeros", {
  expect_true(all(wh_im(M = -1, tau = 1, theta = 2, k = k) == 0))
  expect_true(all(wh_iim(M = -1, tau1 = 0.3, tau0_total = 2, theta = 2, k = k) == 0))
})

test_that("likelihoods are non-negative", {
  expect_true(all(wh_im(M = 3, tau = 2, theta = 1, k = k) >= 0))
  expect_true(all(wh_iim(M = 3, tau1 = 0.5, tau0_total = 4, theta = 1, k = k) >= 0))
  expect_true(all(wh_sc(f = 0.6, tau1 = 0.4, tau0 = 3, theta = 1, k = k) >= 0))
})
