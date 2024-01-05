

test_that(
 desc = 'verbosity prints grow, predict, and importance notes',
 code = {

  skip()

  expect_snapshot(
   fit_verbose <- orsf(pbc, time + status ~.,
                       verbose_progress = TRUE,
                       n_tree = n_tree_test,
                       importance = 'negate')
  )

  expect_snapshot(
   fit_verbose <- orsf(pbc, time + status ~.,
                       verbose_progress = TRUE,
                       n_tree = n_tree_test,
                       importance = 'negate',
                       n_thread = 1)
  )

 }
)

test_that(
 desc = "verbosity is carried by object",
 code = {

  skip()

  fit_verbose <- orsf(penguins, species ~ .,
                      verbose_progress = TRUE,
                      tree_seeds = seeds_standard,
                      n_tree = n_tree_test)

  fit_verbose_2 <- orsf(penguins, bill_length_mm ~ .,
                        verbose_progress = TRUE,
                        tree_seeds = seeds_standard,
                        n_tree = n_tree_test)

  expect_snapshot(
   pd <- orsf_pd_oob(fit_verbose, pred_spec_auto(island))
  )

  expect_snapshot(
   vi <- orsf_vi(fit_verbose_2, importance = 'negate')
  )

  expect_snapshot(
   vs <- orsf_vs(fit_verbose)
  )

 }
)


