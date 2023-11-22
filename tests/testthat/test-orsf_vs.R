
test_that(
 desc = "Survival variable selection filters junk preds",
 code = {

  pbc_with_junk <- pbc

  n_junk_preds <- 5

  junk_names <- paste("junk", seq(n_junk_preds), sep ='_')

  set.seed(329)

  for(i in junk_names)
   pbc_with_junk[[i]] <- rnorm(nrow(pbc))

  fit <- orsf(pbc_with_junk,
              time + status ~ .,
              n_tree = n_tree_test,
              importance = 'anova',
              tree_seeds = seeds_standard)

  fit_var_select <- orsf_vs(fit, n_predictor_min = 3)

  vars_picked <- fit_var_select$predictors_included[[1]]

  expect_false( any(junk_names %in% vars_picked) )

 }
)

test_that(
 desc = "Classification variable selection filters junk preds",
 code = {

  penguins_with_junk <- penguins

  n_junk_preds <- 5

  junk_names <- paste("junk", seq(n_junk_preds), sep ='_')

  set.seed(329)

  for(i in junk_names)
   penguins_with_junk[[i]] <- rnorm(nrow(penguins))

  fit <- orsf(penguins_with_junk,
              species ~ .,
              n_tree = n_tree_test,
              importance = 'permute',
              tree_seeds = seeds_standard)

  fit_var_select <- orsf_vs(fit, n_predictor_min = 3,
                            verbose_progress = TRUE)

  vars_picked <- fit_var_select$predictors_included[[1]]

  expect_false( any(junk_names %in% vars_picked) )

 }
)

test_that(
 desc = "Regression variable selection filters junk preds",
 code = {

  penguins_with_junk <- penguins

  n_junk_preds <- 5

  junk_names <- paste("junk", seq(n_junk_preds), sep ='_')

  set.seed(329)

  for(i in junk_names)
   penguins_with_junk[[i]] <- rnorm(nrow(penguins))

  fit <- orsf(penguins_with_junk,
              species ~ .,
              n_tree = n_tree_test,
              importance = 'permute',
              no_fit = TRUE,
              tree_seeds = seeds_standard)

  fit_var_select <- orsf_vs(fit, n_predictor_min = 3)

  vars_picked <- fit_var_select$predictors_included[[1]]

  expect_false( any(junk_names %in% vars_picked) )

 }
)




