

test_that(
 desc = "fits from orsf_train() are same as orsf()",
 code = {

  skip_on_cran()

  for(i in seq_along(data_list_pbc[c('pbc_standard',
                                     'pbc_status_12',
                                     'pbc_scaled')])){

   for(j in seq_along(controls_surv)){

    fit_untrained <- orsf(data_list_pbc[[i]],
                          formula = time + status ~ . - id,
                          n_tree = n_tree_test,
                          control = controls_surv[[j]],
                          tree_seed = seeds_standard,
                          no_fit = TRUE)

    expect_true(is_empty(fit_untrained$eval_oobag$stat_values))

    expect_false(fit_untrained$trained)

    fit_trained <- orsf_train(fit_untrained)

    expect_equal_leaf_summary(fit_trained, fit_standard_pbc[[j]])

   }

  }

 }

)

