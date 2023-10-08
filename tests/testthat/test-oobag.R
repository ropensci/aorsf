

test_that(
 desc = 'oobag error works w/oobag_eval_every & custom oobag fun works',
 code = {

  fit_custom_oobag <- orsf(pbc,
                           formula = Surv(time, status) ~ .,
                           n_tree = n_tree_test,
                           oobag_eval_every = 1,
                           oobag_fun = oobag_c_survival,
                           tree_seeds = seeds_standard)

  expect_equal_leaf_summary(fit_custom_oobag, fit_standard_pbc$fast)

  expect_equal(
   get_last_oob_stat_value(fit_standard_pbc$fast),
   get_last_oob_stat_value(fit_custom_oobag)
  )

 }
)



