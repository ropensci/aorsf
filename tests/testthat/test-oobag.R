
fit_custom_oobag <- orsf(pbc_orsf,
                         formula = Surv(time, status) ~ . - id,
                         n_tree = 100,
                         oobag_fun = oobag_c_survival,
                         tree_seeds = 1:100)

fit_standard_oobag <- orsf(pbc_orsf,
                           formula = Surv(time, status) ~ . - id,
                           n_tree = 100,
                           tree_seeds = 1:100)

testthat::expect_equal(
 fit_custom_oobag$forest$rows_oobag,
 fit_standard_oobag$forest$rows_oobag
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



