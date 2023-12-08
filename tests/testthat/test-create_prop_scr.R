# Define test data
set.seed(1234)
internal_df <- data.frame(id_col = 1:20, cov1 = rnorm(10, 2), cov2 = rnorm(100, 20))
external_df <- data.frame(id_col = 21:40, cov1 = rnorm(10, 2), cov2 = rnorm(100, 18))
model <- as.formula("~cov1 + cov2")

# Test for valid inputs
test_that("create_prop_scr handles valid inputs", {
  prop_scr_result <- create_prop_scr(internal_df, external_df, id_col, model)
  expect_s3_class(prop_scr_result, "prop_scr")
})

# Test for invalid internal data
test_that("create_prop_scr handles invalid internal data", {
  expect_error(create_prop_scr("abc", external_df, id_col, model))
  expect_error(create_prop_scr(as.matrix(internal_df), external_df, id_col, model))
  expect_error(create_prop_scr(as.tibble(internal_df), external_df, id_col, model))
  expect_error(create_prop_scr(c(5, 10, 15), external_df, id_col, model))
})

# Test for invalid external data
test_that("create_prop_scr handles invalid external data", {
  expect_error(create_prop_scr(internal_df, "abc", id_col, model))
  expect_error(create_prop_scr(internal_df, as.matrix(internal_df), id_col, model))
  expect_error(create_prop_scr(internal_df, as.tibble(internal_df), id_col, model))
  expect_error(create_prop_scr(internal_df, c(5, 10, 15), id_col, model))
})

# Test for invalid model formula
test_that("create_prop_scr handles invalid model formula", {
  expect_error(create_prop_scr(internal_df, external_df, id_col, "formula"))
  expect_error(create_prop_scr(internal_df, external_df, id_col, "cov1+cov2"))
  expect_error(create_prop_scr(internal_df, external_df, id_col, "~cov1+cov2+cov3"))
})

# Test for a missing variable in the data
test_that("create_prop_scr handles missing variables", {
  internal_df_missing_var <- data.frame(id = 1:20, cov1 = rnorm(10, 2))
  expect_error(create_prop_scr(internal_df_missing_var, external_df, id_col, model))
})
