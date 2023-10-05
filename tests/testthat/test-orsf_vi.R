
tree_seeds <- 329

pbc_vi <- pbc_orsf

set.seed(tree_seeds)

pbc_vi$junk <- rnorm(nrow(pbc_orsf))

pbc_vi$junk_cat <- factor(
 sample(letters[1:5], size = nrow(pbc_orsf), replace = TRUE)
)

# simulate a variable with unused factor level
levels(pbc_vi$edema) <- c(levels(pbc_vi$edema), 'empty_lvl')

formula <- Surv(time, status) ~ protime + edema + bili + junk + junk_cat


test_that(
 desc = paste(
  "(1) variable importance is independent from function order",
  "(2) variable importance is independent from n_thread",
  "(3) variable importance is correct"
 ),
 code = {

  for(importance in c('negate', 'permute', 'anova')){

   for(group_factors in c(TRUE, FALSE)){

    fit_with_vi <- orsf(pbc_vi,
                        formula = formula,
                        importance = importance,
                        n_tree = 50,
                        group_factors = group_factors,
                        tree_seeds = tree_seeds)


    vi_during_fit <- orsf_vi(fit_with_vi,
                             group_factors = group_factors)

    if(group_factors){
     expect_true("edema" %in% names(vi_during_fit))
    } else {
     expect_true("edema_0.5" %in% names(vi_during_fit))
     expect_true("edema_1" %in% names(vi_during_fit))
     expect_true(vi_during_fit['edema_empty_lvl'] == 0)
    }


    if(importance != 'anova'){

     fit_no_vi <- orsf(pbc_vi,
                       formula = formula,
                       importance = 'none',
                       n_tree = 50,
                       group_factors = group_factors,
                       tree_seeds = tree_seeds)

     expect_error(orsf_vi(fit_no_vi), regexp = 'no variable importance')

     vi_after_fit <- orsf_vi(fit_no_vi,
                             importance = importance,
                             group_factors = group_factors)

     expect_equal(vi_during_fit, vi_after_fit)

     fit_custom_oobag <- orsf(pbc_vi,
                              formula = formula,
                              importance = importance,
                              n_tree = 50,
                              oobag_fun = oobag_c_risk,
                              group_factors = group_factors,
                              tree_seeds = tree_seeds)

     vi_custom_oobag <- orsf_vi(fit_custom_oobag,
                                group_factors = group_factors)

     # why equal?  oobag_c_risk is a 'custom' eval fun
     # that is equivalent to the eval fun we use by default
     expect_equal(vi_during_fit, vi_custom_oobag)

    }

    fit_threads <- orsf(pbc_vi,
                        formula = formula,
                        importance = importance,
                        n_tree = 50,
                        n_thread = 0,
                        group_factors = group_factors,
                        tree_seeds = tree_seeds)

    vi_threads <- orsf_vi(fit_threads,
                          group_factors = group_factors)

    expect_equal(vi_during_fit, vi_threads)

    good_vars <- c('bili',
                   'protime',
                   if(group_factors) 'edema' else "edema_1")

    bad_vars <- setdiff(names(vi_during_fit), good_vars)

    vi_good_vars <- vi_during_fit[good_vars]
    vi_bad_vars <- vi_during_fit[bad_vars]

    for(j in seq_along(vi_good_vars)){
     expect_true( all(vi_bad_vars < vi_good_vars[j]) )
    }

   }

  }

 }

)

## TODO: move to big output check loop?
# test_that(
#  desc = 'cstat from last run of orsf is reproducible',
#  code = {
#
#   c_target <- last_value(fit$eval_oobag$stat_values)
#   c_estimate <- oobag_c_survival(
#    y_mat = as.matrix(fit$data[, c('time', 'status')]),
#    w_vec = rep(1, nrow(fit$data)),
#    s_vec = fit$pred_oobag
#   )
#
#   expect_equal(c_target, c_estimate)
#
#  }
# )



test_that(
 desc = 'informative errors for custom functions',
 code = {

  expect_error(
   orsf_vi_anova(object = 'nope'),
   regexp = 'inherit'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_bad_name),
   regexp = 'y_mat'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_bad_name_2),
   regexp = 's_vec'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_bad_out),
   regexp = 'length 1'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_bad_out_2),
   regexp = 'type character'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_errors_on_test),
   regexp = 'encountered an error'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_4_args),
   regexp = 'has 4'
  )

 }
)

