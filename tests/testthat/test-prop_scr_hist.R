# Define test data
set.seed(1234)
internal_df <- data.frame(id_col = 1:20, cov1 = rnorm(10, 2), cov2 = rnorm(100, 20))
external_df <- data.frame(id_col = 21:40, cov1 = rnorm(10, 2), cov2 = rnorm(100, 18))
model <- as.formula("~cov1 + cov2")
ps_obj <- calc_prop_scr(internal_df, external_df, id_col, model)

# Test for valid inputs
test_that("prop_scr_hist handles valid inputs", {
  prop_scr_hist_result = prop_scr_hist(ps_obj)
  expect_s3_class(prop_scr_hist_result, "ggplot")
})

# Test for correct output
test_that("prop_scr_hist produces correct output", {
  prop_scr_hist_test <- prop_scr_hist(ps_obj)
  vdiffr::expect_doppelganger("prop-scr-hist-test", prop_scr_hist_test)
})

# Test for invalid inputs
test_that("prop_scr_hist handles invalid inputs", {
  expect_error(prop_scr_hist(6))
  expect_error(prop_scr_hist("k"))
})
