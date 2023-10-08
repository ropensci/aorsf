
test_that(
 desc = "variable selection filters junk preds",
 code = {

  pbc_with_junk <- pbc

  n_junk_preds <- 50

  junk_names <- paste("junk", seq(n_junk_preds), sep ='_')

  set.seed(329)

  for(i in junk_names)
   pbc_with_junk[[i]] <- rnorm(nrow(pbc))

  fit <- orsf(pbc_with_junk,
              time + status ~ .,
              n_tree = 75,
              importance = 'permute',
              tree_seeds = seeds_standard)

  fit_var_select <- orsf_vs(fit, n_predictor_min = 10)

  vars_picked <- fit_var_select$predictors_included[[1]]

  expect_false( any(junk_names %in% vars_picked) )

 }
)




