test_that("Basic Cloud plot",{
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                          external_df = ex_binary_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  basic_plot <- prop_scr_cloud(ps_obj)
  expect_s3_class(basic_plot, "ggplot")
  vdiffr::expect_doppelganger("plot-cloud-test", basic_plot)


})


test_that("Trimmed Cloud plot",{
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                          external_df = ex_binary_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  trimmed_ps_obj <- ps_obj |> trim_ps(low = 0.3)
  trimmed_plot <- prop_scr_cloud(ps_obj, trimmed_ps_obj)

  expect_s3_class(trimmed_plot, "ggplot")
  vdiffr::expect_doppelganger("plot-cloud-trimmed-test", trimmed_plot)



})
