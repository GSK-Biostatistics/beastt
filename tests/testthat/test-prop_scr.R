
test_that("Check is_prop_scr", {
  expect_equal(is_prop_scr(iris), FALSE)
})


test_that("Check trim prop score",{
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                         external_df = ex_binary_df,
                         id_col = subjid,
                         model = ~ cov1 + cov2 + cov3 + cov4)
   trimmed_ps_obj <- trim(ps_obj, low = 0.3, high = 0.7)
   trimmed_df <- trimmed_ps_obj$external_df

   man <- ps_obj$external_df |>
     filter(`___ps___` > 0.3 & `___ps___` < 0.7)
   expect_equal(trimmed_df, man)

  })

test_that("test refitting prop score", {
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                          external_df = ex_binary_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)
  trimmed_ps_obj <- trim(ps_obj, low = 0.3, high = 0.7)
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

  expect_error(rescale(ps_obj))
  expect_error(rescale(ps_obj, n = 75, scale_factor = 1.7))

  pop_recale <- rescale(ps_obj, n = 75)
  # The sum of the new weights should equal the input n
  expect_equal(sum(pop_recale$external_df$`___weight___`),
               75)

  rescale_df <- rescale(ps_obj, scale_factor = 1.2)$external_df

  man <- ps_obj$external_df |>
    mutate(`___weight___` = `___weight___`*1.2)
  expect_equal(rescale_df, man)

})
