# Test for valid inputs
test_that("plot_dist handles valid inputs", {
  plot_dist_result <- plot_dist(dist_normal(10, 2), dist_normal(5, 1))
  expect_s3_class(plot_dist_result, "ggplot")
})

# Test for invalid inputs
test_that("plot_dist handles invalid inputs", {
  expect_error(plot_dist("a"))
  expect_error(plot_dist(5))
})
