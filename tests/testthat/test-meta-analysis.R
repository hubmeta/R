test_that("meta_analysis returns a one-row data frame", {
  out <- meta_analysis(
    correlations = c(.18, .0, .08, .15, .27, .1, .28, .17, .02, .28),
    sample_sizes = c(426, 328, 122, 284, 472, 154, 372, 674, 110, 116),
    reliability_of_x = c(.85, .77, .80, .86, .80, .79, .91, .85, .92, .85),
    reliability_of_y = c(.63, .63, .62, .39, .24, .85, .89, .48, .68, .84),
    significance_levels = c(0.95, 0.80)
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1)
  expect_equal(out$K, 10)
  expect_equal(out$N, 3058)
})
