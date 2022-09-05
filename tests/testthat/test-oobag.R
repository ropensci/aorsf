

fit_custom_oobag <- orsf(pbc_orsf,
                         formula = Surv(time, status) ~ . - id,
                         oobag_fun = oobag_c_harrell,
                         n_tree = 10,
                         tree_seeds = 1:10)

fit_standard_oobag <- orsf(pbc_orsf,
                           formula = Surv(time, status) ~ . - id,
                           n_tree = 10,
                           tree_seeds = 1:10)

test_that(
 desc = 'tree seeds show that a custom oobag fun matches the internal one',
 code = {
  expect_equal(
   fit_standard_oobag$eval_oobag$stat_values,
   fit_custom_oobag$eval_oobag$stat_values
  )
 }
)



