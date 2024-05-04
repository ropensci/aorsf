
test_that(
 desc = "na action of omit works with oobag preds",
 code = {

  mtcars_na <- mtcars

  mtcars_na$vs <- factor(mtcars_na$vs)

  mtcars_na$disp[1] <- NA
  mtcars_na$hp[3] <- NA

  regr_fit <- orsf(
   data = mtcars_na, formula = mpg ~ .,
   n_tree = n_tree_test,
   tree_seeds = seeds_standard,
   na_action = 'omit'
  )

  clsf_fit <- orsf(
   data = mtcars_na, formula = vs ~ .,
   n_tree = n_tree_test,
   tree_seeds = seeds_standard,
   na_action = 'omit'
  )

  surv_fit <- orsf(
   mtcars_na, mpg + vs ~ .,
   n_tree = n_tree_test,
   tree_seeds = seeds_standard,
   na_action = 'omit'
  )

  expect_equal(nrow(na.omit(mtcars_na)), regr_fit$n_obs)
  expect_equal(nrow(na.omit(mtcars_na)), clsf_fit$n_obs)
  expect_equal(nrow(na.omit(mtcars_na)), surv_fit$n_obs)

  clsf_prd_oob <- predict(clsf_fit, oobag = TRUE)
  regr_prd_oob <- predict(regr_fit, oobag = TRUE)
  surv_prd_oob <- predict(surv_fit, oobag = TRUE)

  na_rows <- which(is.na(clsf_fit$pred_oobag))
  expect_true(all(is.na(clsf_prd_oob[na_rows, drop = FALSE])))

  na_rows <- which(is.na(regr_fit$pred_oobag))
  expect_true(all(is.na(regr_prd_oob[na_rows, drop = FALSE])))

  na_rows <- which(is.na(surv_fit$pred_oobag))
  expect_true(all(is.na(surv_prd_oob[na_rows, drop = FALSE])))

 })



test_that(
 desc = "na action of impute works with oobag preds",
 code = {

  mtcars_na <- mtcars

  mtcars_na$vs <- factor(mtcars_na$vs)

  mtcars_na$disp[1] <- NA
  mtcars_na$hp[3] <- NA

  regr_fit <- orsf(
   data = mtcars_na, formula = mpg ~ .,
   n_tree = n_tree_test,
   tree_seeds = seeds_standard,
   na_action = 'impute_meanmode'
  )

  clsf_fit <- orsf(
   data = mtcars_na, formula = vs ~ .,
   n_tree = n_tree_test,
   tree_seeds = seeds_standard,
   na_action = 'impute_meanmode'
  )

  surv_fit <- orsf(
   mtcars_na, mpg + vs ~ .,
   n_tree = n_tree_test,
   tree_seeds = seeds_standard,
   na_action = 'impute_meanmode'
  )

  expect_equal(nrow(mtcars_na), regr_fit$n_obs)
  expect_equal(nrow(mtcars_na), clsf_fit$n_obs)
  expect_equal(nrow(mtcars_na), surv_fit$n_obs)

  clsf_prd_oob <- predict(clsf_fit, oobag = TRUE)
  regr_prd_oob <- predict(regr_fit, oobag = TRUE)
  surv_prd_oob <- predict(surv_fit, oobag = TRUE)

  expect_equal(nrow(mtcars_na), regr_fit$n_obs)
  expect_equal(nrow(mtcars_na), clsf_fit$n_obs)
  expect_equal(nrow(mtcars_na), surv_fit$n_obs)

  clsf_prd_oob <- predict(clsf_fit, oobag = TRUE)
  regr_prd_oob <- predict(regr_fit, oobag = TRUE)
  surv_prd_oob <- predict(surv_fit, oobag = TRUE)

  na_rows <- which(is.na(clsf_fit$pred_oobag))
  expect_true(all(is.na(clsf_prd_oob[na_rows, drop = FALSE])))

  na_rows <- which(is.na(regr_fit$pred_oobag))
  expect_true(all(is.na(regr_prd_oob[na_rows, drop = FALSE])))

  na_rows <- which(is.na(surv_fit$pred_oobag))
  expect_true(all(is.na(surv_prd_oob[na_rows, drop = FALSE])))

 })

