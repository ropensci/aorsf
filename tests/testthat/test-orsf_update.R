

fit <- fit_standard_pbc$fast

fit_new <- orsf_update(fit, formula = . ~ . - edema - spiders,
                       n_tree = n_tree_test + 1)

fit_newer <- orsf_update(fit_new, formula = . ~ . + spiders)

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
