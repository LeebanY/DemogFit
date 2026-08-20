test_that("bsfs_to_s_distribution conserves total block count and matches hand calculation", {
  toy_bsfs <- data.frame(
    count = c(100, 50, 30, 10),
    hetA  = c(0, 1, 0, 0),
    hetB  = c(0, 0, 0, 0),
    hetAB = c(0, 0, 0, 1),
    fixed = c(0, 0, 1, 0)
  )
  s_dist <- bsfs_to_s_distribution(toy_bsfs)

  expect_equal(sum(s_dist$count), sum(toy_bsfs$count))
  # hand-computed expected result (see package development notes)
  expect_equal(s_dist$count[s_dist$k == 0], 130)
  expect_equal(s_dist$count[s_dist$k == 1], 60)
})
