
# survival tests ----------------------------------------------------------

test_that(
 desc = paste(
  "Survival forest:",
  "(1) variable importance is independent from function order",
  "(2) variable importance is independent from n_thread",
  "(3) variable importance is correct",
  collapse = '\n'
 ),
 code = {

  skip_on_cran()

  pbc_vi <- pbc_orsf

  pbc_vi$junk <- rnorm(nrow(pbc_orsf))

  pbc_vi$junk_cat <- factor(
   sample(letters[1:5], size = nrow(pbc_orsf), replace = TRUE)
  )

  # simulate a variable with unused factor level
  levels(pbc_vi$edema) <- c(levels(pbc_vi$edema), 'empty_lvl')

  formula <- Surv(time, status) ~ protime + edema + bili + junk + junk_cat

  for(importance in c('negate', 'permute', 'anova')){

   for(group_factors in c(TRUE, FALSE)){

    fit_with_vi <- orsf(pbc_vi,
                        formula = formula,
                        importance = importance,
                        n_tree = n_tree_test,
                        group_factors = group_factors,
                        tree_seeds = seeds_standard)


    vi_during_fit <- orsf_vi(fit_with_vi, group_factors = group_factors)

    wrapper_fun <- switch(
     importance,
     'anova' = orsf_vi_anova,
     'permute' = orsf_vi_permute,
     'negate' = orsf_vi_negate
    )

    expect_equal(
     vi_during_fit,
     wrapper_fun(fit_with_vi, group_factors = group_factors)
    )


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
                       n_tree = n_tree_test,
                       group_factors = group_factors,
                       tree_seeds = seeds_standard)

     expect_error(orsf_vi(fit_no_vi), regexp = 'no variable importance')

     vi_after_fit <- orsf_vi(fit_no_vi,
                             importance = importance,
                             group_factors = group_factors)

     fit_vi_custom <- orsf(pbc_vi,
                           formula = formula,
                           n_tree = n_tree_test,
                           oobag_fun = oobag_c_risk,
                           importance = importance,
                           tree_seeds = seeds_standard)

     vi_custom_during_fit <- orsf_vi(fit_vi_custom,
                                     group_factors = group_factors)

     vi_custom_after_fit <- orsf_vi(fit_no_vi,
                                    importance = importance,
                                    group_factors = group_factors,
                                    oobag_fun = oobag_c_risk)


     expect_equal(vi_during_fit, vi_after_fit)
     expect_equal(vi_custom_after_fit, vi_after_fit)
     expect_equal(vi_custom_during_fit, vi_after_fit)

     fit_custom_oobag <- orsf(pbc_vi,
                              formula = formula,
                              importance = importance,
                              n_tree = n_tree_test,
                              oobag_fun = oobag_c_risk,
                              group_factors = group_factors,
                              tree_seeds = seeds_standard)

     vi_custom_oobag <- orsf_vi(fit_custom_oobag,
                                group_factors = group_factors)

     # why equal?  oobag_c_risk is a 'custom' eval fun
     # that is equivalent to the eval fun we use by default
     expect_equal(vi_during_fit, vi_custom_oobag)

    }

    fit_threads <- orsf(pbc_vi,
                        formula = formula,
                        importance = importance,
                        n_tree = n_tree_test,
                        n_thread = 0,
                        group_factors = group_factors,
                        tree_seeds = seeds_standard)

    vi_threads <- orsf_vi(fit_threads,
                          group_factors = group_factors)

    expect_equal(vi_during_fit, vi_threads)

    good_vars <- c('bili', 'protime')

    bad_vars <- c('junk',
                  if(group_factors) 'junk_cat'
                  else paste("junk_cat", levels(pbc_vi$junk_cat)[-1], sep = '_'))

    vi_good_vars <- vi_during_fit[good_vars]
    vi_bad_vars <- vi_during_fit[bad_vars]

    for(j in seq_along(vi_good_vars)){
     expect_true( mean(vi_bad_vars < vi_good_vars[j]) > 1/2 )


     if(mean(vi_bad_vars < vi_good_vars[j]) <= 1/2){

      print("WFT")
      print(vi_good_vars)
      print(vi_bad_vars)

     }

    }

   }

  }

  fit_no_data <- orsf(pbc_vi,
                      formula = time + status ~ .,
                      n_tree = n_tree_test,
                      attach_data = FALSE)

  expect_error(orsf_vi(fit_no_data), regexp = 'training data')

 }

)

# classification tests -----------------------------------------------------

test_that(
 desc = paste(
  "Classification forest:",
  "(1) variable importance is independent from function order",
  "(2) variable importance is independent from n_thread",
  "(3) variable importance is correct",
  collapse = '\n'
 ),
 code = {

  penguins_vi <- penguins_orsf

  penguins_vi$junk <- rnorm(nrow(penguins_orsf))

  penguins_vi$junk_cat <- factor(
   sample(letters[1:5], size = nrow(penguins_orsf), replace = TRUE)
  )

  # simulate a variable with unused factor level
  levels(penguins_vi$island) <- c(levels(penguins_vi$island), 'empty_lvl')

  formula <- species ~ .

  for(importance in c('negate', 'permute', 'anova')){

   for(group_factors in c(TRUE, FALSE)){

    fit_with_vi <- orsf(penguins_vi,
                        formula = formula,
                        importance = importance,
                        n_tree = n_tree_test,
                        group_factors = group_factors,
                        tree_seeds = seeds_standard)

    vi_during_fit <- orsf_vi(fit_with_vi,
                             group_factors = group_factors)

    wrapper_fun <- switch(
     importance,
     'anova' = orsf_vi_anova,
     'permute' = orsf_vi_permute,
     'negate' = orsf_vi_negate
    )

    expect_equal(
     vi_during_fit,
     wrapper_fun(fit_with_vi, group_factors = group_factors)
    )

    if(group_factors){
     expect_true("island" %in% names(vi_during_fit))
    } else {
     expect_true("island_Dream" %in% names(vi_during_fit))
     expect_true("island_Torgersen" %in% names(vi_during_fit))
     expect_true("island_empty_lvl" %in% names(vi_during_fit))
     expect_true(vi_during_fit['island_empty_lvl'] == 0)
    }

    if(importance != 'anova'){

     fit_no_vi <- orsf(penguins_vi,
                       formula = formula,
                       importance = 'none',
                       n_tree = n_tree_test,
                       group_factors = group_factors,
                       tree_seeds = seeds_standard)

     expect_error(orsf_vi(fit_no_vi), regexp = 'no variable importance')

     vi_after_fit <- orsf_vi(fit_no_vi,
                             importance = importance,
                             group_factors = group_factors)

     fit_vi_custom <- orsf(penguins_vi,
                           formula = formula,
                           n_tree = n_tree_test,
                           oobag_fun = oobag_brier_clsf,
                           importance = importance,
                           tree_seeds = seeds_standard)

     vi_custom_during_fit <- orsf_vi(fit_vi_custom,
                                     group_factors = group_factors)

     vi_custom_after_fit <- orsf_vi(fit_no_vi,
                                    importance = importance,
                                    group_factors = group_factors,
                                    oobag_fun = oobag_brier_clsf)


     expect_equal(vi_during_fit, vi_after_fit)
     expect_equal(vi_custom_during_fit, vi_custom_after_fit)

    }

    fit_threads <- orsf(penguins_vi,
                        formula = formula,
                        importance = importance,
                        n_tree = n_tree_test,
                        n_thread = 0,
                        group_factors = group_factors,
                        tree_seeds = seeds_standard)

    vi_threads <- orsf_vi(fit_threads,
                          group_factors = group_factors)

    expect_equal(vi_during_fit, vi_threads)

    good_vars <- c('bill_length_mm',
                   'flipper_length_mm',
                   'body_mass_g')

    bad_vars <- c('junk',
                  if(group_factors) 'junk_cat'
                  else paste("junk_cat", levels(penguins_vi$junk_cat)[-1], sep = '_'))

    vi_good_vars <- vi_during_fit[good_vars]
    vi_bad_vars <- vi_during_fit[bad_vars]

    for(j in seq_along(vi_good_vars)){
     expect_true( mean(vi_bad_vars < vi_good_vars[j]) > 1/2 )
    }

   }

  }

  fit_no_vi <- orsf(penguins_vi,
                    formula = formula,
                    n_tree = n_tree_test,
                    importance = 'none')
  expect_error(orsf_vi_anova(fit_no_vi), regexp = 'ANOVA')
  expect_error(orsf_vi(fit_no_vi, importance = 'anova'), regexp = 'ANOVA')

 }
)

# regression tests --------------------------------------------------------

test_that(
 desc = paste(
  "Classification forest:",
  "(1) variable importance is independent from function order",
  "(2) variable importance is independent from n_thread",
  "(3) variable importance is correct",
  collapse = '\n'
 ),
 code = {

  skip_on_cran()

  penguins_vi <- penguins_orsf

  penguins_vi$junk <- rnorm(nrow(penguins_orsf))

  penguins_vi$junk_cat <- factor(
   sample(letters[1:5], size = nrow(penguins_orsf), replace = TRUE)
  )

  # simulate a variable with unused factor level
  levels(penguins_vi$island) <- c(levels(penguins_vi$island), 'empty_lvl')

  # still using penguin data, but switching the outcome
  formula <- bill_length_mm ~ .


  for(importance in c('negate', 'permute', 'anova')){

   for(group_factors in c(TRUE, FALSE)){

    fit_with_vi <- orsf(penguins_vi,
                        formula = formula,
                        importance = importance,
                        n_tree = n_tree_test,
                        group_factors = group_factors,
                        tree_seeds = seeds_standard)

    vi_during_fit <- orsf_vi(fit_with_vi, group_factors = group_factors)

    wrapper_fun <- switch(
     importance,
     'anova' = orsf_vi_anova,
     'permute' = orsf_vi_permute,
     'negate' = orsf_vi_negate
    )

    expect_equal(
     vi_during_fit,
     wrapper_fun(fit_with_vi, group_factors = group_factors)
    )

    if(group_factors){
     expect_true("island" %in% names(vi_during_fit))
    } else {
     expect_true("island_Dream" %in% names(vi_during_fit))
     expect_true("island_Torgersen" %in% names(vi_during_fit))
     expect_true("island_empty_lvl" %in% names(vi_during_fit))
     expect_true(vi_during_fit['island_empty_lvl'] == 0)
    }

    if(importance != 'anova'){

     fit_no_vi <- orsf(penguins_vi,
                       formula = formula,
                       importance = 'none',
                       n_tree = n_tree_test,
                       group_factors = group_factors,
                       tree_seeds = seeds_standard)

     expect_error(orsf_vi(fit_no_vi), regexp = 'no variable importance')

     vi_after_fit <- orsf_vi(fit_no_vi,
                             importance = importance,
                             group_factors = group_factors)

     # oobag_fun looks like a typo here, but it is not a typo.
     # oobag_brier_clsf is equivalent to regression R-squared

     fit_vi_custom <- orsf(penguins_vi,
                           formula = formula,
                           n_tree = n_tree_test,
                           oobag_fun = oobag_brier_clsf,
                           importance = importance,
                           tree_seeds = seeds_standard)

     vi_custom_during_fit <- orsf_vi(fit_vi_custom,
                                     group_factors = group_factors)

     vi_custom_after_fit <- orsf_vi(fit_no_vi,
                                    importance = importance,
                                    group_factors = group_factors,
                                    oobag_fun = oobag_brier_clsf)


     expect_equal(vi_during_fit, vi_after_fit)
     expect_equal(vi_custom_during_fit, vi_custom_after_fit)
     expect_equal(vi_after_fit, vi_custom_after_fit)

    }

    fit_threads <- orsf(penguins_vi,
                        formula = formula,
                        importance = importance,
                        n_tree = n_tree_test,
                        n_thread = 0,
                        group_factors = group_factors,
                        tree_seeds = seeds_standard)

    vi_threads <- orsf_vi(fit_threads,
                          group_factors = group_factors)

    expect_equal(vi_during_fit, vi_threads)

    good_vars <- c('flipper_length_mm',
                   if(group_factors) 'species'
                   else c('species_Chinstrap', 'species_Gentoo'),
                   'body_mass_g')

    bad_vars <- c('junk',
                  if(group_factors) 'junk_cat'
                  else paste("junk_cat", levels(penguins_vi$junk_cat)[-1], sep = '_'))

    vi_good_vars <- vi_during_fit[good_vars]
    vi_bad_vars <- vi_during_fit[bad_vars]

    for(j in seq_along(vi_good_vars)){
     expect_true( mean(vi_bad_vars < vi_good_vars[j]) > 1/2 )
    }

   }

  }

 }

)



# General tests -----------------------------------------------------------

test_that(
 desc = 'informative errors for custom functions',
 code = {

  skip_on_cran()

  fit_no_vi <- orsf_update(fit_standard_pbc$fast, importance = 'none')

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
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_bad_name_3),
   regexp = 'w_vec'
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

test_that(
 desc = 'calling convenience vi functions does not change forest data',
 code = {

  importance_init <- fit_standard_pbc$fast$importance
  importance_type_init <- fit_standard_pbc$fast$importance_type

  importance_permute <- orsf_vi_permute(fit_standard_pbc$fast)

  expect_equal(importance_init, fit_standard_pbc$fast$importance)
  expect_equal(importance_type_init, fit_standard_pbc$fast$importance_type)

  importance_negate <- orsf_vi_negate(fit_standard_pbc$fast)

  expect_equal(importance_init, fit_standard_pbc$fast$importance)
  expect_equal(importance_type_init, fit_standard_pbc$fast$importance_type)

 }
)

