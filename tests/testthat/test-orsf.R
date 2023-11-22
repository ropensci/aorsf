
test_that(
 desc = 'non-formula inputs are vetted',
 code = {

  # correct formula
  f <- time + status ~ .

  expect_error(orsf(pbc, f, n_tree = 0), "should be >= 1")
  expect_error(orsf(pbc, f, n_split = "3"), "should have type")
  expect_error(orsf(pbc, f, mtry = 5000), 'should be <=')
  expect_error(orsf(pbc, f, leaf_min_events = 5000), 'should be <=')
  expect_error(orsf(pbc, f, leaf_min_obs = 5000), 'should be <=')
  expect_error(orsf(pbc, f, attachData = TRUE), 'attach_data?')
  expect_error(orsf(pbc, f, Control = 0), 'control?')
  expect_error(orsf(pbc, f, tree_seeds = c(1,2,3)), 'number of trees')
  expect_error(orsf(pbc, f, sample_fraction = 1, oobag_pred_type = 'risk'),
               'no samples are out-of-bag')
  expect_error(orsf(pbc, f, split_rule = 'cstat', split_min_stat = 1),
               'should be < 1')

  pbc_orsf$date_var <- Sys.Date()
  expect_error(orsf(pbc_orsf, f), 'unsupported type')
  pbc_orsf$date_var <- NULL

 }
)

test_that(
 desc = 'outcome type can be guessed',
 code = {


  fit_regr <- orsf(penguins, bill_length_mm ~ ., no_fit = TRUE)
  fit_clsf <- orsf(penguins, species ~ ., no_fit = TRUE)
  fit_surv <- orsf(pbc, time + status ~ ., no_fit = TRUE)

  expect_s3_class(fit_regr, "ObliqueForestRegression")
  expect_s3_class(fit_clsf, "ObliqueForestClassification")
  expect_s3_class(fit_surv, "ObliqueForestSurvival")

 }
)

test_that(
 desc = 'potential user-errors with outcome types are caught',
 code = {

  expect_error(
   orsf(penguins, species ~., control = orsf_control_regression()),
   "it is a factor"
  )

  expect_error(
   orsf(penguins, bill_length_mm ~., control = orsf_control_classification()),
   "please convert bill_length_mm to a factor"
  )

 }
)


test_that(
 desc = 'target_df too high is caught',
 code = {

  cntrl <- orsf_control_survival(method = 'net', target_df = 10)
  expect_error(orsf(pbc, time + status ~ ., control = cntrl), 'should be <=')

 }
)

test_that(
 desc = 'orsf runs the same with data.table vs. data.frame',
 code = {

  fit_dt <- orsf(as.data.table(pbc),
                 formula = time + status ~ .,
                 n_tree = n_tree_test,
                 control = controls_surv$fast,
                 tree_seed = seeds_standard)

  expect_equal_leaf_summary(fit_dt, fit_standard_pbc$fast)

  fit_dt <- orsf(as.data.table(penguins),
                 formula = species ~ .,
                 n_tree = n_tree_test,
                 control = controls_clsf$fast,
                 tree_seed = seeds_standard)

  expect_equal_leaf_summary(fit_dt, fit_standard_penguin_species$fast)

  fit_dt <- orsf(as.data.table(penguins),
                 formula = bill_length_mm ~ .,
                 n_tree = n_tree_test,
                 control = controls_regr$fast,
                 tree_seed = seeds_standard)

  expect_equal_leaf_summary(fit_dt, fit_standard_penguin_bills$fast)

 }
)

test_that(
 desc = "orsf runs with lists and recipes",
 code = {

  pbc_list <- as.list(pbc_orsf)
  pbc_list_bad <- pbc_list
  pbc_list_bad$trt <- pbc_list_bad$trt[1:3]
  pbc_list_bad$age <- pbc_list_bad$age[1:5]

  # skip() # I don't want to list recipes in suggests
  #
  # recipe <- recipes::recipe(pbc_orsf, formula = time + status ~ .) %>%
  #  recipes::step_rm(id)
  #
  # recipe_prepped <- recipes::prep(recipe)
  #
  # fit_recipe <- orsf(recipe_prepped, Surv(time, status) ~ .,
  #                    n_tree = n_tree_test,
  #                    tree_seeds = seeds_standard)
  #
  # expect_equal_leaf_summary(fit_recipe, fit_standard_pbc$fast)

  fit_list <- orsf(pbc_list,
                   Surv(time, status) ~ . - id,
                   n_tree = n_tree_test,
                   tree_seeds = seeds_standard)

  expect_equal_leaf_summary(fit_list, fit_standard_pbc$fast)

  expect_error(
   orsf(pbc_list_bad, Surv(time, status) ~ .),
   regexp = 'unable to cast data'
  )

 }
)

test_that(
 desc = "blank and non-standard names trigger an error",
 code = {

  pbc_temp <- pbc
  pbc_temp$x1 <- rnorm(nrow(pbc_temp))
  pbc_temp$x2 <- rnorm(nrow(pbc_temp))

  names(pbc_temp)[names(pbc_temp)=='x1'] <- ""
  names(pbc_temp)[names(pbc_temp)=='x2'] <- " "

  expect_error(
   orsf(data = pbc_temp, Surv(time, status) ~ . - id), regex = 'Blank'
  )


  pbc_temp <- pbc
  pbc_temp$x1 <- rnorm(nrow(pbc_temp))
  pbc_temp$x2 <- rnorm(nrow(pbc_temp))

  names(pbc_temp)[names(pbc_temp)=='x1'] <- "@"
  names(pbc_temp)[names(pbc_temp)=='x2'] <- "#"

  expect_error(
   orsf(data = pbc_temp, Surv(time, status) ~ . - id), regex = 'Non\\-standard'
  )

 }
)

test_that(
 desc = 'if oobag time is unspecified, pred horizon = median(time)',
 code = {

  fit_1 <- orsf(data = pbc_orsf,
                formula = time + status ~ . - id,
                n_tree = 1)

  fit_2 <- orsf(data = pbc_orsf,
                formula = time + status ~ . - id,
                n_tree = 1,
                oobag_pred_type = 'none')

  expect_equal(fit_1$pred_horizon, median(pbc_orsf$time))
  expect_equal(fit_1$pred_horizon, fit_2$pred_horizon)

 }
)


test_that(
 desc = 'list columns are not allowed',
 code = {

  pbc_temp <- pbc_orsf
  pbc_temp$list_col <- list(list(a=1))

  expect_error(
   orsf(pbc_temp, time + status ~ . - id),
   regexp = '<list_col>'
  )
 }
)


test_that(
 desc = "algorithm grows more accurate with higher number of iterations",
 code = {

  n_tree <- n_tree_test * 5
  eval_every <- max(round(n_tree/5), 1)

  fit <- orsf(pbc,
              formula = Surv(time, status) ~ .,
              n_tree = n_tree,
              leaf_min_obs = 50,
              tree_seeds = seeds_standard,
              oobag_eval_every = eval_every)

  expect_lt(fit$eval_oobag$stat_values[1],
            last_value(fit$eval_oobag$stat_values))

  fit <- orsf(penguins,
              formula = species ~ .,
              n_tree = n_tree,
              leaf_min_obs = 50,
              tree_seeds = seeds_standard,
              oobag_eval_every = eval_every)

  expect_lt(fit$eval_oobag$stat_values[1],
            last_value(fit$eval_oobag$stat_values))

  fit <- orsf(penguins,
              formula = bill_length_mm ~ .,
              leaf_min_obs = 50,
              n_tree = n_tree, # just needs a bit extra
              tree_seeds = seeds_standard,
              oobag_eval_every = eval_every)

  expect_lt(fit$eval_oobag$stat_values[1],
            last_value(fit$eval_oobag$stat_values))

 }
)


test_that(
 desc = 'Empty training data throw an error',
 code = {

  expect_error(
   orsf(pbc_orsf[c(), ], Surv(time, status) ~ . - id),
   regexp = 'training data are empty'
  )

  expect_error(
   orsf(pbc_orsf[, c()], Surv(time, status) ~ . - id),
   regexp = 'training data are empty'
  )

 }
)

test_that(
 desc = "Data with all-`NA` fields or columns are rejected",
 code = {

  pbc_temp <- pbc
  pbc_temp[, 'bili'] <- NA_real_

  expect_error(orsf(pbc_temp, time + status ~ . - id,
                    na_action = 'omit'),
               'complete data')

  expect_error(orsf(pbc_temp, time + status ~ . - id,
                    na_action = 'impute_meanmode'),
               'column bili has no observed values')

 }
)

test_that(
 desc = "data with missing values are rejected when na_action is fail",
 code = {

  pbc_temp <- pbc
  pbc_temp[1, 'bili'] <- NA_real_

  expect_error(orsf(pbc_temp, time + status ~ . - id),
               'missing values')

 }
)

test_that(
 desc = 'missing data are dropped when na_action is omit',
 code = {

  pbc_temp <- pbc
  pbc_temp[1, 'bili'] <- NA_real_

  fit_omit <- orsf(pbc_temp, time + status ~ .-id,  na_action = 'omit')

  expect_equal(fit_omit$n_obs, nrow(stats::na.omit(pbc_temp)))

 }
)


test_that(
 desc = 'robust to threading, outcome formats, scaling, and noising',
 code = {

  fits_surv <- lapply(data_list_pbc[-1],  function(data){
   orsf(data,
        formula = time + status ~ .,
        n_thread = 2,
        n_tree = n_tree_test,
        tree_seeds = seeds_standard)
  })

  expect_equal_leaf_summary(fits_surv$pbc_status_12,
                            fit_standard_pbc$fast)

  expect_equal_oobag_eval(fits_surv$pbc_scaled,
                          fit_standard_pbc$fast,
                          tolerance = .01)

  expect_equal_oobag_eval(fits_surv$pbc_noised,
                          fit_standard_pbc$fast,
                          tolerance = .01)

  fits_clsf <- lapply(data_list_penguins[-1],  function(data){
   orsf(data,
        formula = species ~ .,
        n_thread = 2,
        n_tree = n_tree_test,
        tree_seeds = seeds_standard)
  })

  expect_equal(fits_clsf$penguins_binary$n_class, 2)
  expect_equal(fits_clsf$penguins_scaled$n_class, 3)

  expect_equal(ncol(fits_clsf$penguins_binary$pred_oobag), 2)

  expect_equal_oobag_eval(fits_clsf$penguins_scaled,
                          fit_standard_penguin_species$fast,
                          tolerance = .01)

  expect_equal_oobag_eval(fits_clsf$penguins_noised,
                          fit_standard_penguin_species$fast,
                          tolerance = .01)

  fits_regr <- lapply(data_list_penguins[-1],  function(data){
   orsf(data,
        formula = bill_length_mm ~ .,
        n_thread = 2,
        n_tree = n_tree_test,
        tree_seeds = seeds_standard)
  })

  expect_equal_oobag_eval(fits_regr$penguins_scaled,
                          fit_standard_penguin_bills$fast,
                          tolerance = .01)

  expect_equal_oobag_eval(fits_regr$penguins_scaled,
                          fit_standard_penguin_bills$fast,
                          tolerance = .01)

 }
)


test_that(
 desc = 'oob error correct for user-specified function',
 code = {

  fit <- orsf(data = pbc,
              formula = time + status ~ . -id,
              n_tree = n_tree_test,
              oobag_fun = oobag_c_risk,
              tree_seeds = seeds_standard)

  expect_equal_oobag_eval(fit, fit_standard_pbc$fast)

  # can also reproduce it from the oobag predictions
  expect_equal(
   oobag_c_risk(
    y_mat = as.matrix(pbc_orsf[,c("time", "status")]),
    w_vec = rep(1, nrow(pbc_orsf)),
    s_vec = fit$pred_oobag
   ),
   as.numeric(fit$eval_oobag$stat_values)
  )

  # skip() # don't want to suggest yardstick or Hmisc
  #
  # oobag_rsq_eval <- function(y_mat, w_vec, s_vec){
  #
  #  yardstick::rsq_trad_vec(truth = as.numeric(y_mat),
  #                          estimate = as.numeric(s_vec),
  #                          case_weights = as.numeric(w_vec))
  # }
  #
  # fit <- orsf(data = mtcars,
  #             formula = mpg ~ .,
  #             n_tree = n_tree_test,
  #             oobag_fun = oobag_rsq_eval,
  #             tree_seeds = seeds_standard)
  #
  # expect_equal(
  #  fit$eval_oobag$stat_values[1,1],
  #  yardstick::rsq_trad_vec(truth = as.numeric(mtcars$mpg),
  #                          estimate = as.numeric(fit$pred_oobag),
  #                          case_weights = rep(1, nrow(mtcars)))
  # )
  #
  # oobag_cstat_clsf <- function(y_mat, w_vec, s_vec){
  #
  #  y_vec = as.numeric(y_mat)
  #  cstat <- Hmisc::somers2(x = s_vec,
  #                          y = y_vec,
  #                          weights = w_vec)['C']
  #  cstat
  #
  # }
  #
  # fit <- orsf(data = penguins,
  #             formula = species ~ .,
  #             n_tree = n_tree_test,
  #             oobag_fun = oobag_cstat_clsf,
  #             tree_seeds = seeds_standard)
  #
  # expect_equal_oobag_eval(fit, fit_standard_penguin_species$fast)


 }
)


test_that(
 desc = 'orsf_time_to_train is reasonable at approximating time to train',
 code = {

  # testing the seed behavior when no_fit is TRUE. You should get the same
  # forest whether you train with orsf() or with orsf_train().

  object <- orsf(pbc, Surv(time, status) ~ .,
                 n_tree = n_tree_test,
                 tree_seeds = 1,
                 no_fit = TRUE,
                 importance = 'none')

  time_estimated <- orsf_time_to_train(object, n_tree_subset = 1)

  time_true_start <- Sys.time()
  fit_orsf_3 <- orsf_train(object)
  time_true_stop <- Sys.time()

  time_true <- time_true_stop - time_true_start

  diff <- abs(as.numeric(time_true - time_estimated))

  # estimated time is within 5 seconds of true time.
  expect_lt(diff, 5)

 }
)

test_that(
 desc = 'orsf_fit objects can be saved and loaded with saveRDS and readRDS',
 code = {

  skip_on_cran()

  fil <- tempfile("fit_orsf", fileext = ".rds")

  ## save a single object to file
  saveRDS(fit_standard_pbc$fast, fil)
  ## restore it under a different name
  fit <- readRDS(fil)

  p1 <- predict(fit_standard_pbc$fast, new_data = pbc_test)
  p2 <- predict(fit, new_data = pbc_test)

  expect_equal(p1, p2)

 }
)

test_that(
 desc = 'weights do not make trees grow more than intended',
 code = {

  fit_unwtd <- orsf(pbc, time + status ~ .,
                    n_tree = n_tree_test,
                    tree_seeds = seeds_standard)

  fit_wtd <- orsf(pbc,
                  time + status ~ .,
                  weights = rep(2, nrow(pbc_orsf)),
                  n_tree = n_tree_test,
                  tree_seeds = seeds_standard)

  # using weights should not inadvertently make trees deeper.
  expect_equal(fit_wtd$get_mean_leaves_per_tree(),
               fit_unwtd$get_mean_leaves_per_tree(),
               tolerance = 1/2)

 }
)



test_that(
 desc = 'user-supplied beta functions are vetted',
 code = {

  f_bad_1 <- function(a_node, y_node, w_node){ 1 }
  f_bad_2 <- function(x_node, a_node, w_node){ 1 }
  f_bad_3 <- function(x_node, y_node, a_node){ 1 }
  f_bad_4 <- function(x_node, y_node){ 1 }

  f_bad_5 <- function(x_node, y_node, w_node) {
   stop("an expected error occurred")
  }

  f_bad_6 <- function(x_node, y_node, w_node){
   return(matrix(0, ncol = 2, nrow = ncol(x_node)))
  }

  f_bad_7 <- function(x_node, y_node, w_node){
   return(matrix(0, ncol = 1, nrow = 2))
  }

  f_bad_8 <- function(x_node, y_node, w_node) {runif(n = ncol(x_node))}

  expect_error(
   orsf(pbc, time + status ~ .,
        control = orsf_control_survival(method = f_bad_1)),
   'x_node'
  )
  expect_error(
   orsf(pbc, time + status ~ .,
        control = orsf_control_survival(method = f_bad_2)),
   'y_node'
  )
  expect_error(
   orsf(pbc, time + status ~ .,
        control = orsf_control_survival(method = f_bad_3)),
   'w_node'
  )
  expect_error(
   orsf(pbc, time + status ~ .,
        control = orsf_control_survival(method = f_bad_4)),
   'should have 3'
  )
  expect_error(
   orsf(pbc, time + status ~ .,
        control = orsf_control_survival(method = f_bad_5)),
   'encountered an error'
  )
  expect_error(
   orsf(pbc, time + status ~ .,
        control = orsf_control_survival(method = f_bad_6)),
   'with 1 column'
  )
  expect_error(
   orsf(pbc, time + status ~ .,
        control = orsf_control_survival(method = f_bad_7)),
   'with 1 row for each'
  )
  expect_error(
   orsf(pbc, time + status ~ .,
        control = orsf_control_survival(method = f_bad_8)),
   'matrix output'
  )

 }
)

test_that(
 desc = "user supplied beta functions are applied correctly",
 code = {

  fit_pca = orsf(pbc,
                 Surv(time, status) ~ .,
                 tree_seeds = seeds_standard,
                 control = orsf_control_survival(method = f_pca),
                 n_tree = n_tree_test)

  expect_gt(fit_pca$eval_oobag$stat_values, .765)

 }
)

test_that(
 desc = 'oblique survival forests run as intended for valid inputs',
 code = {

  # just takes forever.
  skip_on_cran()

  inputs <- expand.grid(
   data_format = c('plain'),
   n_tree = 1,
   n_split = 1,
   n_retry = 0,
   mtry = 3,
   sample_with_replacement = c(TRUE, FALSE),
   leaf_min_events = 5,
   leaf_min_obs = c(10),
   split_rule = c("logrank", "cstat"),
   split_min_events = 5,
   split_min_obs = 15,
   oobag_pred_type = c('none', 'risk', 'mort'),
   oobag_pred_horizon = c(1,2,3),
   orsf_control = c('cph', 'net', 'custom'),
   stringsAsFactors = FALSE
  )

  for(i in seq(nrow(inputs))){

   data_fun <- switch(
    as.character(inputs$data_format[i]),
    'plain' = function(x) x,
    'tibble' = tibble::as_tibble,
    'data.table' = as.data.table
   )

   pred_horizon <- switch(inputs$oobag_pred_horizon[i],
                          '1' = 1000,
                          '2' = c(1000, 2000),
                          '3' = c(1000, 2000, 3000))

   control <- switch(inputs$orsf_control[i],
                     'cph' = orsf_control_survival(method = 'glm'),
                     'net' = orsf_control_survival(method = 'net'),
                     'custom' = orsf_control_survival(method = f_pca))

   if(inputs$sample_with_replacement[i]){
    sample_fraction <- 0.632
   } else {
    sample_fraction <- runif(n = 1, min = .25, max = .75)
   }

   fit <- orsf(data = data_fun(pbc_orsf),
               formula = time + status ~ . - id,
               control = control,
               sample_with_replacement = inputs$sample_with_replacement[i],
               sample_fraction = sample_fraction,
               n_tree = inputs$n_tree[i],
               n_split = inputs$n_split[i],
               n_retry = inputs$n_retry[i],
               mtry = inputs$mtry[i],
               leaf_min_events = inputs$leaf_min_events[i],
               leaf_min_obs = inputs$leaf_min_obs[i],
               split_rule = inputs$split_rule[i],
               split_min_events = inputs$split_min_events[i],
               split_min_obs = inputs$split_min_obs[i],
               oobag_pred_type = inputs$oobag_pred_type[i],
               oobag_pred_horizon = pred_horizon)

   expect_s3_class(fit, class = 'ObliqueForestSurvival')

   # data are not unintentionally modified by reference,
   expect_identical(data_fun(pbc_orsf), fit$data)


   expect_no_missing(fit$forest)
   expect_no_missing(fit$importance)
   expect_no_missing(fit$pred_horizon)

   expect_length(fit$forest$rows_oobag,   n = fit$n_tree)
   expect_length(fit$forest$cutpoint,     n = fit$n_tree)
   expect_length(fit$forest$child_left,   n = fit$n_tree)
   expect_length(fit$forest$coef_indices, n = fit$n_tree)
   expect_length(fit$forest$coef_values,  n = fit$n_tree)
   expect_length(fit$forest$leaf_summary, n = fit$n_tree)

   if(!inputs$sample_with_replacement[i]){
    expect_equal(
     1 - length(fit$forest$rows_oobag[[1]]) / fit$n_obs,
     sample_fraction,
     tolerance = 0.025
    )
   }

   if(inputs$oobag_pred_type[i] != 'none'){

    if(inputs$oobag_pred_type[i] %in% c("chf","surv","risk")){

     expect_length(fit$eval_oobag$stat_values, length(pred_horizon))

    } else if(inputs$oobag_pred_type[i] == 'mort'){

     expect_length(fit$eval_oobag$stat_values, 1)

    }


    expect_equal(nrow(fit$pred_oobag), fit$n_obs)

    # these lengths should match for n_tree=1
    # b/c only the oobag rows of the first tree
    # will get a prediction value. Note that the
    # vectors themselves aren't equal b/c rows_oobag
    # corresponds to the sorted version of the data.
    expect_equal(
     length(which(complete.cases(fit$pred_oobag))),
     length(fit$forest$rows_oobag[[1]])
    )

    oobag_preds <- na.omit(fit$pred_oobag)

    expect_true(all(oobag_preds >= 0))

    if(inputs$oobag_pred_type[i] %in% c("risk", "surv")){
     expect_true(all(oobag_preds <= 1))
    }

   } else {
    expect_equal(dim(fit$eval_oobag$stat_values), c(0, 0))
   }

  }

 }
)

test_that(
 desc = 'oblique classification forests run as intended for valid inputs',
 code = {

  # just takes forever.
  skip_on_cran()

  inputs <- expand.grid(
   data_format = c('plain'),
   n_tree = 1,
   n_split = 1,
   n_retry = 0,
   mtry = 3,
   sample_with_replacement = c(TRUE, FALSE),
   leaf_min_obs = 10,
   split_rule = c("gini", "cstat"),
   split_min_obs = 15,
   oobag_pred_type = c('none', 'prob'),
   orsf_control = c('glm', 'net', 'custom'),
   stringsAsFactors = FALSE
  )

  for(i in seq(nrow(inputs))){

   data_fun <- switch(
    as.character(inputs$data_format[i]),
    'plain' = function(x) x,
    'tibble' = tibble::as_tibble,
    'data.table' = as.data.table
   )

   control <- switch(inputs$orsf_control[i],
                     'glm' = orsf_control_classification(method = 'glm'),
                     'net' = orsf_control_classification(method = 'net'),
                     'custom' = orsf_control_classification(method = f_pca))

   if(inputs$sample_with_replacement[i]){
    sample_fraction <- 0.632
   } else {
    sample_fraction <- runif(n = 1, min = .25, max = .75)
   }

   fit <- orsf(data = data_fun(penguins_orsf),
               formula = species ~ .,
               control = control,
               sample_with_replacement = inputs$sample_with_replacement[i],
               sample_fraction = sample_fraction,
               n_tree = inputs$n_tree[i],
               n_split = inputs$n_split[i],
               n_retry = inputs$n_retry[i],
               mtry = inputs$mtry[i],
               leaf_min_events = inputs$leaf_min_events[i],
               leaf_min_obs = inputs$leaf_min_obs[i],
               split_rule = inputs$split_rule[i],
               split_min_events = inputs$split_min_events[i],
               split_min_obs = inputs$split_min_obs[i],
               oobag_pred_type = inputs$oobag_pred_type[i])

   expect_s3_class(fit, class = 'ObliqueForestClassification')

   # data are not unintentionally modified by reference,
   expect_identical(data_fun(penguins_orsf), fit$data)


   expect_no_missing(fit$forest)
   expect_no_missing(fit$importance)

   expect_length(fit$forest$rows_oobag,   n = fit$n_tree)
   expect_length(fit$forest$cutpoint,     n = fit$n_tree)
   expect_length(fit$forest$child_left,   n = fit$n_tree)
   expect_length(fit$forest$coef_indices, n = fit$n_tree)
   expect_length(fit$forest$coef_values,  n = fit$n_tree)
   expect_length(fit$forest$leaf_summary, n = fit$n_tree)

   if(!inputs$sample_with_replacement[i]){
    expect_equal(
     1 - length(fit$forest$rows_oobag[[1]]) / fit$n_obs,
     sample_fraction,
     tolerance = 0.025
    )
   }

   if(inputs$oobag_pred_type[i] != 'none'){

    expect_length(fit$eval_oobag$stat_values, 1)

    expect_equal(nrow(fit$pred_oobag), fit$n_obs)

    # these lengths should match for n_tree=1
    # b/c only the oobag rows of the first tree
    # will get a prediction value. Note that the
    # vectors themselves aren't equal b/c rows_oobag
    # corresponds to the sorted version of the data.
    expect_equal(
     length(which(complete.cases(fit$pred_oobag))),
     length(fit$forest$rows_oobag[[1]])
    )

    oobag_preds <- na.omit(fit$pred_oobag)

    expect_true(all(apply(oobag_preds, 1, sum) == 1))
    expect_true(all(oobag_preds >= 0))
    expect_true(all(oobag_preds <= 1))

   } else {

    expect_equal(dim(fit$eval_oobag$stat_values), c(0, 0))

   }

  }

 }
)

test_that(
 desc = 'oblique regression forests run as intended for valid inputs',
 code = {

  # just takes forever.
  skip_on_cran()

  inputs <- expand.grid(
   data_format = c('plain'),
   n_tree = 1,
   n_split = 1,
   n_retry = 0,
   mtry = 3,
   sample_with_replacement = c(TRUE, FALSE),
   leaf_min_obs = 3,
   split_rule = c("variance"),
   split_min_obs = 6,
   oobag_pred_type = c('none', 'mean'),
   orsf_control = c('glm', 'net', 'custom'),
   stringsAsFactors = FALSE
  )

  for(i in seq(nrow(inputs))){

   data_fun <- switch(
    as.character(inputs$data_format[i]),
    'plain' = function(x) x,
    'tibble' = tibble::as_tibble,
    'data.table' = as.data.table
   )

   control <- switch(inputs$orsf_control[i],
                     'glm' = orsf_control_regression(method = 'glm'),
                     'net' = orsf_control_regression(method = 'net'),
                     'custom' = orsf_control_regression(method = f_pca))

   if(inputs$sample_with_replacement[i]){
    sample_fraction <- 0.632
   } else {
    sample_fraction <- runif(n = 1, min = .25, max = .75)
   }

   fit <- orsf(data = data_fun(penguins),
               formula = bill_length_mm ~ .,
               control = control,
               sample_with_replacement = inputs$sample_with_replacement[i],
               sample_fraction = sample_fraction,
               n_tree = inputs$n_tree[i],
               n_split = inputs$n_split[i],
               n_retry = inputs$n_retry[i],
               mtry = inputs$mtry[i],
               leaf_min_events = inputs$leaf_min_events[i],
               leaf_min_obs = inputs$leaf_min_obs[i],
               split_rule = inputs$split_rule[i],
               split_min_events = inputs$split_min_events[i],
               split_min_obs = inputs$split_min_obs[i],
               oobag_pred_type = inputs$oobag_pred_type[i])

   expect_s3_class(fit, class = 'ObliqueForestRegression')

   # data are not unintentionally modified by reference,
   expect_identical(data_fun(penguins), fit$data)


   expect_no_missing(fit$forest)
   expect_no_missing(fit$importance)

   expect_length(fit$forest$rows_oobag,   n = fit$n_tree)
   expect_length(fit$forest$cutpoint,     n = fit$n_tree)
   expect_length(fit$forest$child_left,   n = fit$n_tree)
   expect_length(fit$forest$coef_indices, n = fit$n_tree)
   expect_length(fit$forest$coef_values,  n = fit$n_tree)
   expect_length(fit$forest$leaf_summary, n = fit$n_tree)

   if(!inputs$sample_with_replacement[i]){
    expect_equal(
     1 - length(fit$forest$rows_oobag[[1]]) / fit$n_obs,
     sample_fraction,
     # bigger tolerance b/c sample size is small
     tolerance = 0.075
    )
   }

   if(inputs$oobag_pred_type[i] != 'none'){

    expect_length(fit$eval_oobag$stat_values, 1)

    expect_equal(nrow(fit$pred_oobag), fit$n_obs)

    # these lengths should match for n_tree=1
    # b/c only the oobag rows of the first tree
    # will get a prediction value. Note that the
    # vectors themselves aren't equal b/c rows_oobag
    # corresponds to the sorted version of the data.
    expect_equal(
     length(which(complete.cases(fit$pred_oobag))),
     length(fit$forest$rows_oobag[[1]])
    )


   } else {

    expect_equal(dim(fit$eval_oobag$stat_values), c(0, 0))

   }

  }

 }
)
