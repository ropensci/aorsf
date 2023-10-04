

test_that(
 desc = "untrained forest acts the same as a trained one",
 code = {

  fit_untrained <- orsf(pbc_orsf,
                        formula = time + status ~ . - id,
                        n_tree = 10,
                        tree_seed = seeds_standard,
                        no_fit = TRUE)

  expect_true(is_empty(fit_untrained$eval_oobag$stat_values))


  expect_equal( attr(fit_untrained, 'trained'), FALSE )

  fit_trained <- orsf_train(fit_untrained)

  expect_equal(fit_trained$forest$leaf_summary,
               fit_standard$forest$leaf_summary)


 }
)

