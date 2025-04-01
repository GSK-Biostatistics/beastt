
test_that("Check is_prop_scr", {
  expect_equal(is_prop_scr(iris), FALSE)
})


test_that("Check trim prop score",{
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                         external_df = ex_binary_df,
                         id_col = subjid,
                         model = ~ cov1 + cov2 + cov3 + cov4)

  # Direct
   trimmed_ps_obj <- trim_ps(ps_obj, low = 0.3, high = 0.7)
   trimmed_df <- trimmed_ps_obj$external_df

   man <- ps_obj$external_df |>
     filter(`___ps___` > 0.3 & `___ps___` < 0.7)
   expect_equal(trimmed_df, man)

   # Test with only low boundary
   low_only <- trim_ps(ps_obj, low = 0.3)
   expect_true(all(low_only$external_df$`___ps___` >= 0.3))

   # Test with only high boundary
   high_only <- trim_ps(ps_obj, high = 0.7)
   expect_true(all(high_only$external_df$`___ps___` <= 0.7))

   # Manual calculation for comparison
   man_low_only <- ps_obj$external_df |>
     filter(`___ps___` >= 0.3)
   expect_equal(low_only$external_df, man_low_only)

   man_high_only <- ps_obj$external_df |>
     filter(`___ps___` <= 0.7)
   expect_equal(high_only$external_df, man_high_only)

   # Quantile
   trimmed_df_high <- trim_ps(ps_obj, high = 0.75, quantile = TRUE)$external_df
   trimmed_df_low <- trim_ps(ps_obj, low = 0.25,  quantile = TRUE)$external_df

   ps_vals <- ps_obj$external_df |>
     pull(`___ps___` )
   low_cv <- quantile(ps_vals, 0.25)
   high_cv <-  quantile(ps_vals, 0.75)

   man_low <-  ps_obj$external_df |>
     filter(`___ps___` >= low_cv)
   expect_equal(trimmed_df_low, man_low)

   man_high <-  ps_obj$external_df |>
     filter(`___ps___` <= high_cv)
   expect_equal(trimmed_df_high, man_high)


   # Errors
   expect_error(trim_ps(ps_obj, low = -0.3))
   expect_error(trim_ps(ps_obj, high = 1.2))


  })


test_that("trim preserves prop_scr object structure", {
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                          external_df = ex_binary_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  trimmed_ps_obj <- trim_ps(ps_obj, low = 0.2, high = 0.8)

  # Check that required properties exist in trimmed object
  expect_true(is_prop_scr(trimmed_ps_obj))
  expect_equal(names(ps_obj), names(trimmed_ps_obj))
  expect_equal(class(ps_obj), class(trimmed_ps_obj))

  # Check model formula is preserved
  expect_equal(ps_obj$model_formula, trimmed_ps_obj$model_formula)

  # Check that id column is preserved
  expect_equal(ps_obj$id_col, trimmed_ps_obj$id_col)
})

test_that("trim with NULL parameters returns unchanged dataset", {
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                          external_df = ex_binary_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  # Trim with NULL parameters should return the original object
  null_trim <- trim_ps(ps_obj, low = NULL, high = NULL)
  expect_equal(ps_obj$external_df, null_trim$external_df)
})

test_that("trim correctly handles quantile edge cases", {
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                          external_df = ex_binary_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  # Test with quantile = 0
  q0_trim <- trim_ps(ps_obj, low = 0, quantile = TRUE)
  expect_equal(nrow(q0_trim$external_df), nrow(ps_obj$external_df))

  # Test with quantile = 1
  q1_trim <- trim_ps(ps_obj, high = 1, quantile = TRUE)
  expect_equal(nrow(q1_trim$external_df), nrow(ps_obj$external_df))

  # Test with very small quantile range
  narrow_trim <- trim_ps(ps_obj, low = 0.49, high = 0.51, quantile = TRUE)
  ps_vals <- ps_obj$external_df |>
    pull(`___ps___`)
  low_val <- quantile(ps_vals, 0.49)
  high_val <- quantile(ps_vals, 0.51)

  expect_true(all(narrow_trim$external_df$`___ps___` >= low_val))
  expect_true(all(narrow_trim$external_df$`___ps___` <= high_val))
})


test_that("test refitting prop score", {
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                          external_df = ex_binary_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)
  trimmed_ps_obj <- trim_ps(ps_obj, low = 0.3, high = 0.7)
  # Manual calc
  dat <- bind_rows(trimmed_ps_obj$external_df,
            trimmed_ps_obj$internal_df)



  # Calculating the absolute standardized mean difference
  asmd_adj <- bal.tab(select(dat,cov1, cov2, cov3, cov4 ), # df of covariates (internal and external)
                      treat = dat$`___internal___`,   # internal indicator
                      binary = "std",         # use standardized version of mean differences for binary covariates
                      continuous = "std",     # use standardized version of mean differences for continuous covariates
                      s.d.denom = "pooled",   # calculation of the denominator of SMD
                      weights = dat$`___weight___`,
                      abs = TRUE)$Balance

  asmd_unadj <- bal.tab(select(dat,cov1, cov2, cov3, cov4 ), # df of covariates (internal and external)
                        treat = dat$`___internal___`,   # internal indicator
                        binary = "std",         # use standardized version of mean differences for binary covariates
                        continuous = "std",     # use standardized version of mean differences for continuous covariates
                        s.d.denom = "pooled",   # calculation of the denominator of SMD
                        abs = TRUE)$Balance

  asmd_man <- tibble(
    covariate = rownames(asmd_adj),
    diff_unadj = asmd_unadj[,2],
    diff_adj = asmd_adj[,3],
  )

  expect_equal(trimmed_ps_obj$abs_std_mean_diff, asmd_man)

  })

test_that("Check rescale prop score",{
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                          external_df = ex_binary_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  expect_error(rescale_ps(ps_obj))
  expect_error(rescale_ps(ps_obj, n = 75, scale_factor = 1.7))

  pop_recale <- rescale_ps(ps_obj, n = 75)
  # The sum of the new weights should equal the input n
  expect_equal(sum(pop_recale$external_df$`___weight___`),
               75)

  rescale_df <- rescale_ps(ps_obj, scale_factor = 1.2)$external_df

  man <- ps_obj$external_df |>
    mutate(`___weight___` = `___weight___`*1.2)
  expect_equal(rescale_df, man)

})
