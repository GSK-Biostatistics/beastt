# Define test data
set.seed(1234)
internal_df <- data.frame(id_col = 1:20, cov1 = rnorm(10, 2), cov2 = rnorm(100, 20))
external_df <- data.frame(id_col = 21:40, cov1 = rnorm(10, 2), cov2 = rnorm(100, 18))
model <- as.formula("~cov1 + cov2")
ps_obj <- calc_prop_scr(internal_df, external_df, id_col, model)

# Test for valid inputs
test_that("prop_scr_dens handles valid inputs", {
  prop_scr_dens_result = prop_scr_dens(ps_obj, variable="ipw")
  expect_s3_class(prop_scr_dens_result, "ggplot")
})

# Test for invalid inputs
test_that("prop_scr_dens handles invalid inputs", {
  expect_error(prop_scr_dens(25))
  expect_error(prop_scr_dens("cat"))
})
