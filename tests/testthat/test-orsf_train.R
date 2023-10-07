

test_that(
 desc = "fits from orsf_train() are same as orsf()",
 code = {

  for(i in seq_along(data_list_pbc[c('pbc_standard',
                                     'pbc_status_12',
                                     'pbc_scaled')])){

   for(j in seq_along(controls)){

    fit_untrained <- orsf(data_list_pbc[[i]],
                          formula = time + status ~ . - id,
                          n_tree = n_tree_test,
                          control = controls[[j]],
                          tree_seed = seeds_standard,
                          no_fit = TRUE)

    expect_true(is_empty(fit_untrained$eval_oobag$stat_values))

    expect_equal( attr(fit_untrained, 'trained'), FALSE )

    fit_trained <- orsf_train(fit_untrained)

    expect_equal_leaf_summary(fit_trained, fit_standard_pbc[[j]])

   }

  }

 }

)

