
sink("orsf_output_1.txt")
fit_custom_oobag <- orsf(pbc_orsf,
                         formula = Surv(time, status) ~ . - id,
                         n_tree = 100,
                         tree_seeds = 1:100,
                         verbose_progress = 0)
sink()
sink("orsf_output_2.txt")
fit_standard_oobag <- orsf(pbc_orsf,
                           formula = Surv(time, status) ~ . - id,
                           n_tree = 100,
                           tree_seeds = 1:100,
                           verbose_progress = 0)
sink()

testthat::expect_equal(
 fit_custom_oobag$forest$rows_oobag,
 fit_standard_oobag$forest$rows_oobag
)

testthat::expect_equal(
 fit_custom_oobag$forest$cutpoint,
 fit_standard_oobag$forest$cutpoint
)


test_that(
 desc = 'tree seeds show that a custom oobag fun matches the internal one',
 code = {
  expect_equal(
   fit_standard_oobag$eval_oobag$stat_values,
   fit_custom_oobag$eval_oobag$stat_values
  )
 }
)



