
# classification ----

# double weights
w_dbl <- runif(nrow(penguins_orsf))
# integer weights
w_int <- sample(1:3, nrow(penguins_orsf), replace = TRUE)

fit_unwt <- fit_standard_penguin_species$fast

fit_w_dbl <- orsf(species ~ ., data = penguins_orsf,
                  tree_seeds = seeds_standard,
                  n_tree = n_tree_test,
                  weights = w_dbl)

fit_w_int <- orsf(species ~ ., data = penguins_orsf,
                  tree_seeds = seeds_standard,
                  n_tree = n_tree_test,
                  weights = w_int)

pd_unwt  <- orsf_pd_oob(fit_unwt, pred_spec_auto(bill_length_mm))
pd_w_dbl <- orsf_pd_oob(fit_w_dbl, pred_spec_auto(bill_length_mm))
pd_w_int <- orsf_pd_oob(fit_w_int, pred_spec_auto(bill_length_mm))

vi_unwt  <- orsf_vi(fit_unwt)
vi_w_dbl <- orsf_vi(fit_w_dbl)
vi_w_int <- orsf_vi(fit_w_int)

test_that(
 desc = "weights do not impact randomness",
 code = {
  expect_equal(fit_w_int$forest$rows_oobag,
               fit_w_dbl$forest$rows_oobag)
  expect_equal(fit_unwt$forest$rows_oobag,
               fit_w_dbl$forest$rows_oobag)
 }
)

test_that(
 desc = 'weights are used in pd',
 code = {
  expect_true(all(pd_unwt$mean != pd_w_int$mean))
  expect_true(all(pd_w_dbl$mean != pd_w_int$mean))
 }
)

test_that(
 desc = 'weights are used in vi',
 code = {
  expect_true(any(vi_unwt != vi_w_int))
  expect_true(any(vi_w_dbl != vi_w_int))
 }
)

test_that(
 desc = 'dropping weights gives same forest as never having them',
 code = {

  fit_unwt_update <- orsf_update(fit_w_int,
                                 weights = rep(1, nrow(penguins_orsf)))

  expect_equal_leaf_summary(fit_unwt, fit_unwt_update)

 }
)
