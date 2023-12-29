

fit <- fit_standard_pbc$fast

fit_new <- orsf_update(fit, formula = . ~ . - edema - spiders,
                       n_tree = n_tree_test + 1)

fit_newer <- orsf_update(fit_new, formula = . ~ . + spiders)

fit_untrained <- orsf_update(fit, no_fit = TRUE)

fit_trained <- orsf_update(fit_untrained, no_fit = FALSE)

test_that(
 desc = "train and untrain update trained status and parameters",
 code = {


  expect_false(fit_untrained$trained)
  expect_null(fit_untrained$get_means())
  expect_equal(fit_untrained$get_mean_leaves_per_tree(), 0)

  expect_true(fit_trained$trained)
  expect_equal(fit_trained$get_means(), fit$get_means())
  expect_equal(fit_trained$get_mean_leaves_per_tree(),
               fit$get_mean_leaves_per_tree())
 }
)

test_that(
 desc = 'variables can be removed and added with updated formula',
 code = {
  expect_false("edema" %in% all.vars(fit_new$formula))
  expect_false("spiders" %in% all.vars(fit_new$formula))
  expect_true("spiders" %in% all.vars(fit_newer$formula))
 }
)

pbc_new <- select_cols(pbc_test, all.vars(fit_new$formula))

test_that(
 desc = "error handlers are updated with updated fits",
 code = {
  expect_error(predict(fit, new_data = pbc_new),
               regexp = 'columns not contained')
  expect_type(predict(fit_new, new_data = pbc_new), 'double')
 }
)

test_that(
 desc = 'updated forest has same tree seeds on overlap',
 code = {
  expect_equal(fit$tree_seeds, fit_new$tree_seeds[seq(n_tree_test)])
 }
)

test_that(
 desc = "dot checker catches the dots",
 code = {
  expect_error(orsf_update(fit, ntree = 10),
               regexp = 'did you mean n_tree?')
 }
)

test_that(
 desc = "fit from scratch is equivalent to fit from update",
 code = {

  fit_vet_from_scratch <- orsf(survival::veteran,
                               time + status ~ .,
                               n_tree = n_tree_test,
                               tree_seeds = seeds_standard)

  fit_vet_from_update <- orsf_update(fit, data = survival::veteran)

  expect_equal_leaf_summary(fit_vet_from_scratch, fit_vet_from_update)

 }
)

test_that(
 desc = "setting a default null field to null reverts that field to the default value",
 code = {

  fit_control <- fit_standard_penguin_bills$custom

  fit_control_null <- orsf_update(fit_control, control = NULL)

  expect_equal_leaf_summary(fit_standard_penguin_bills$fast,
                            fit_control_null)

  fit_split <- fit_standard_pbc$fast

  fit_split_c <- orsf_update(fit_split, split_rule = 'cstat')

  expect_true(fit_split_c$split_rule == 'cstat')

  fit_revert <- orsf_update(fit_split_c, split_rule = NULL)

  expect_true(fit_revert$split_rule == 'logrank')

  expect_equal_leaf_summary(fit_revert, fit_split)

 }
)


