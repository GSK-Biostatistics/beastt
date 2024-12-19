# Define test data
set.seed(1234)
internal_df <- data.frame(id_col = 1:20, cov1 = rnorm(10, 2), cov2 = rnorm(100, 20))
external_df <- data.frame(id_col = 21:40, cov1 = rnorm(10, 2), cov2 = rnorm(100, 18))
model <- as.formula("~cov1 + cov2")
ps_obj <- calc_prop_scr(internal_df, external_df, id_col, model)

# Test for valid inputs
test_that("prop_scr_love handles valid inputs", {
  prop_scr_love_result <- prop_scr_love(ps_obj)
  expect_s3_class(prop_scr_love_result, "ggplot")
})

# Test for correct output
test_that("prop_scr_love produces correct output", {
  prop_scr_love_test <- prop_scr_love(ps_obj)
  vdiffr::expect_doppelganger("prop-scr-love-test", prop_scr_love_test)
})

# Test for invalid inputs
test_that("prop_scr_love handles invalid inputs", {
  expect_error(prop_scr_love("love"))
  expect_error(plot_scr_love(4))
})
