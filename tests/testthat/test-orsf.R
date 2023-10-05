

f <- time + status ~ . - id

test_that(
 desc = 'non-formula inputs are vetted',
 code = {

  expect_error(orsf(pbc_orsf, f, n_tree = 0), "should be >= 1")
  expect_error(orsf(pbc_orsf, f, n_split = "3"), "should have type")
  expect_error(orsf(pbc_orsf, f, mtry = 5000), 'should be <=')
  expect_error(orsf(pbc_orsf, f, leaf_min_events = 5000), 'should be <=')
  expect_error(orsf(pbc_orsf, f, leaf_min_obs = 5000), 'should be <=')
  expect_error(orsf(pbc_orsf, f, attachData = TRUE), 'attach_data?')
  expect_error(orsf(pbc_orsf, f, Control = 0), 'control?')

  pbc_orsf$date_var <- Sys.Date()
  expect_error(orsf(pbc_orsf, f), 'unsupported type')
  pbc_orsf$date_var <- NULL

 }
)

test_that(
 desc = 'target_df too high is caught',
 code = {

  cntrl <- orsf_control_net(df_target = 10)
  expect_error(orsf(pbc_orsf, formula = f, control = cntrl), 'must be <= mtry')

 }
)

test_that(
 desc = 'orsf runs with data.table and with net control',
 code = {

  expect_s3_class(orsf(as.data.table(pbc_orsf), f, n_tree = 1), 'orsf_fit')

  expect_s3_class(orsf(as.data.table(pbc_orsf), f,
                       control = orsf_control_net(),
                       n_tree = 3), 'orsf_fit')
 }
)


#' @srrstats {G5.8b, G5.8b} *Data of unsupported types trigger an error*

test_that(
 desc = "blank and non-standard names trigger an error",
 code = {

  pbc_temp <- pbc_orsf
  pbc_temp$x1 <- rnorm(nrow(pbc_temp))
  pbc_temp$x2 <- rnorm(nrow(pbc_temp))

  names(pbc_temp)[names(pbc_temp)=='x1'] <- ""
  names(pbc_temp)[names(pbc_temp)=='x2'] <- " "

  expect_error(
   orsf(data = pbc_temp, Surv(time, status) ~ . - id), regex = 'Blank'
  )


  pbc_temp <- pbc_orsf
  pbc_temp$x1 <- rnorm(nrow(pbc_temp))
  pbc_temp$x2 <- rnorm(nrow(pbc_temp))

  names(pbc_temp)[names(pbc_temp)=='x1'] <- "@"
  names(pbc_temp)[names(pbc_temp)=='x2'] <- "#"

  expect_error(
   orsf(data = pbc_temp, Surv(time, status) ~ . - id), regex = 'Non\\-standard'
  )

 }
)


#' @srrstats {G2.11} *testing allowance and accounting for units class*

test_that(
 'orsf tracks meta data for units class variables',
 code = {

  suppressMessages(library(units))
  pbc_units <- pbc_orsf


  units(pbc_units$time) <- 'days'
  units(pbc_units$age) <- 'years'
  units(pbc_units$bili) <- 'mg/dl'

  fit_units <- orsf(pbc_units, Surv(time, status) ~ . - id, n_tree=1)

  expect_equal(
   get_unit_info(fit_units),
   list(
    time = list(
     numerator = "d",
     denominator = character(0),
     label = "d"
    ),
    age = list(
     numerator = "years",
     denominator = character(0),
     label = "years"
    ),
    bili = list(
     numerator = "mg",
     denominator = "dl",
     label = "mg/dl"
    )
   )
  )
 }

)

data_fit <- copy(pbc_orsf)

fit_with_vi <- orsf(data = data_fit,
                    formula = Surv(time, status) ~ . - id,
                    importance = 'negate',
                    n_tree = 50)

test_that("data are not unintentionally modified by reference",
          code = {expect_identical(data_fit, pbc_orsf)})


fit_no_vi <- orsf(data = pbc_orsf,
                  formula = Surv(time, status) ~ . - id,
                  importance = 'none',
                  n_tree = 50)


#' @srrstats {G5.3} *Explicit test expected to return objects containing no missing (`NA`) or undefined (`NaN`, `Inf`) values are explicitly tested.*

test_that(
 "output contains no missing values",
 code = {

  miss_check_no_vi <- sapply(fit_no_vi, no_miss_list)
  miss_check_with_vi <- sapply(fit_with_vi, no_miss_list)

  for(i in seq_along(miss_check_no_vi)){
   if(!is_empty(miss_check_no_vi[[i]])){

    if(is.matrix(miss_check_no_vi[[i]])){
     miss_check_no_vi[[i]] <- unlist(miss_check_no_vi[[i]])
    }
    expect_true(sum(miss_check_no_vi[[i]]) == 0)

   }

  }

  for(i in seq_along(miss_check_with_vi)){
   if(!is_empty(miss_check_with_vi[[i]])){

    if(is.matrix(miss_check_with_vi[[i]])){
     miss_check_with_vi[[i]] <- unlist(miss_check_with_vi[[i]])
    }
    expect_true(sum(miss_check_with_vi[[i]]) == 0)

   }
  }


 }
)


#' @srrstats {G5.7} **Algorithm performance tests** *test that implementation performs as expected as properties of data change. These tests shows that as data size increases, fit time increases. Conversely, fit time decreases as convergence thresholds increase. Also, fit time decreases as the maximum iterations decrease.*

# I'm making the difference in data size very big because I don't want this
# test to fail on some operating systems.
pbc_small <- pbc_orsf[1:50, ]

pbc_large <- rbind(pbc_orsf, pbc_orsf, pbc_orsf, pbc_orsf, pbc_orsf)
pbc_large <- rbind(pbc_large, pbc_large, pbc_large, pbc_large, pbc_large)

test_that(
 desc = "algorithm runs slower as data size increases",
 code = {
  time_small <- system.time(orsf(pbc_small,
                                 Surv(time, status) ~ . -id,
                                 n_tree=50))

  time_large <- system.time(orsf(pbc_large,
                                 Surv(time, status) ~ . -id,
                                 n_tree=50))

  expect_true(time_small['elapsed'] < time_large['elapsed'])
 }
)


test_that(
 desc = "algorithm runs faster with lower convergence tolerance",
 code = {

  time_small <- system.time(
   orsf(pbc_orsf,
        control = orsf_control_fast(),
        Surv(time, status) ~ . -id,
        n_tree = 500)
  )

  time_large <- system.time(
   orsf(pbc_orsf,
        control = orsf_control_cph(iter_max = 50, eps = 1e-10),
        Surv(time, status) ~ . -id,
        n_tree = 500)
  )

  expect_true(time_small['elapsed'] < time_large['elapsed'])

 }
)

test_that(
 desc = "algorithm runs faster with lower number of iterations",
 code = {

  time_small <- system.time(
   orsf(pbc_orsf,
        Surv(time, status) ~ . -id,
        n_tree = 5)
  )

  time_large <- system.time(
   orsf(pbc_orsf,
        Surv(time, status) ~ . -id,
        n_tree = 1000) # big difference prevents unneeded failure
  )

  expect_true(time_small['elapsed'] < time_large['elapsed'])

 }
)


#' @srrstats {ML7.11} *OOB C-statistic is monitored by this test. As the number of trees in the forest increases, the C-statistic should also increase*

test_that(
 desc = "algorithm grows more accurate with higher number of iterations",
 code = {

  fit <- orsf(pbc_orsf,
              formula = Surv(time, status) ~ . -id,
              n_tree = 50,
              oobag_eval_every = 5)

  expect_lt(fit$eval_oobag$stat_values[1],
            last_value(fit$eval_oobag$stat_values))

 }
)


#' @srrstats {G5.8, G5.8a} **Edge condition tests** *Zero-length data produce expected behaviour*

test_that(
 desc = 'Boundary case: empty training data throw an error',
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

#' @srrstats {G5.8c, G5.8c} *Data with all-`NA` fields or columns are rejected*

pbc_temp <- pbc_orsf
pbc_temp[, 'bili'] <- NA_real_

test_that(
 desc = "Data with all-`NA` fields or columns are rejected",
 code = {
  expect_error(orsf(pbc_temp, time + status ~ . - id,
                    na_action = 'omit'),
               'column bili has no observed values')

  expect_error(orsf(pbc_temp, time + status ~ . - id,
                    na_action = 'impute_meanmode'),
               'column bili has no observed values')

 }
)

pbc_temp$bili[1:10] <- 12

test_that(
 desc = "data with missing values are rejected when na_action is fail",
 code = {

  expect_error(orsf(pbc_temp, time + status ~ . - id),
               'missing values')


 }
)

pbc_temp <- copy(pbc_orsf)
pbc_temp[1:10, 'bili'] <- NA_real_
pbc_temp_orig <- copy(pbc_temp)

test_that(
 desc = 'missing data are dropped when na_action is omit',
 code = {

  fit_omit <- orsf(pbc_temp, time + status ~ .-id,  na_action = 'omit')
  expect_equal(nrow(fit_omit$data),
               nrow(stats::na.omit(pbc_temp)))

 }
)

test_that(
 desc = 'missing data are imputed when na_action is impute_meanmode',
 code = {

  fit_impute <- orsf(pbc_temp,
                     time + status ~ .,
                     na_action = 'impute_meanmode')

  expect_equal(fit_impute$data$bili[1:10],
               rep(mean(pbc_temp$bili, na.rm=TRUE), 10))

  expect_equal(fit_impute$data$bili[-c(1:10)],
               pbc_temp$bili[-c(1:10)])

 }
)


test_that("data are not unintentionally modified by reference when imputed",
          code = {expect_identical(pbc_temp, pbc_temp_orig)})

#' @srrstats {G5.9} **Noise susceptibility tests**
#' @srrstats {G5.9a} *Adding trivial noise to data does not meaningfully change results*
#' @srrstats {G5.9b} *Running under different random seeds gives identifal results*

pbc_noise <- pbc_orsf
pbc_scale <- pbc_orsf

vars <- c('bili', 'chol', 'albumin', 'copper', 'alk.phos', 'ast')

set.seed(730)

for(i in vars){
 pbc_noise[[i]] <- add_noise(pbc_noise[[i]])
 pbc_scale[[i]] <- change_scale(pbc_scale[[i]])
}


fit_orsf <-
 orsf(pbc_orsf, Surv(time, status) ~ . - id,
      n_thread = 1,
      n_tree = 500,
      tree_seeds = 1:500)

fit_orsf_2 <-
 orsf(pbc_orsf, Surv(time, status) ~ . - id,
      n_thread = 5,
      n_tree = 500,
      tree_seeds = 1:500)

fit_orsf_noise <-
 orsf(pbc_noise, Surv(time, status) ~ . - id,
      n_tree = 500,
      tree_seeds = 1:500)

fit_orsf_scale <-
 orsf(pbc_scale, Surv(time, status) ~ . - id,
      n_tree = 500,
      tree_seeds = 1:500)

#' @srrstats {ML7.1} *Demonstrate effect of numeric scaling of input data.*
test_that(
 desc = 'outputs are robust to multi-threading, scaling, and noising',
 code = {

  expect_lt(
   abs(
    fit_orsf$eval_oobag$stat_values -
     fit_orsf_scale$eval_oobag$stat_values
   ),
   0.01
  )

  expect_lt(
   abs(
    fit_orsf$eval_oobag$stat_values -
     fit_orsf_2$eval_oobag$stat_values
   ),
   0.01
  )

  expect_lt(
   abs(
    fit_orsf$eval_oobag$stat_values -
     fit_orsf_noise$eval_oobag$stat_values
   ),
   0.01
  )


  expect_lt(
   max(abs(fit_orsf$pred_oobag - fit_orsf_scale$pred_oobag)),
   0.1
  )

  expect_lt(
   max(abs(fit_orsf$pred_oobag - fit_orsf_2$pred_oobag)),
   0.1
  )

  expect_lt(
   max(abs(fit_orsf$pred_oobag - fit_orsf_noise$pred_oobag)),
   0.1
  )

  expect_lt(
   mean(abs(fit_orsf$importance - fit_orsf_noise$importance)),
   0.1
  )

  expect_equal(fit_orsf$forest,
               fit_orsf_2$forest)

  expect_equal(fit_orsf$importance,
               fit_orsf_2$importance)

  expect_equal(fit_orsf$forest$rows_oobag,
               fit_orsf_noise$forest$rows_oobag)

  expect_equal(fit_orsf$forest$rows_oobag,
               fit_orsf_scale$forest$rows_oobag)

  expect_equal(fit_orsf$forest$leaf_summary,
               fit_orsf_scale$forest$leaf_summary)

 }
)

test_that(
 desc = 'results are identical if a forest is fitted under the same random seed',
 code = {

  object <- orsf(pbc_orsf, Surv(time, status) ~ . - id,
                 n_tree = 500,
                 tree_seeds = 1:500,
                 no_fit = TRUE)
  fit_orsf_3 <- orsf_train(object)

  expect_equal(fit_orsf$forest,
               fit_orsf_3$forest)

  attr_orsf <- attributes(fit_orsf)
  attr_orsf_3 <- attributes(fit_orsf_3)

  for(i in names(attr_orsf)){

   if( !(i %in% c('f_beta', 'f_oobag_eval')) ){

    expect_equal(attr_orsf[[i]], attr_orsf_3[[i]])

   }

  }

 }

)

test_that(
 desc = 'oob rows identical with same tree seeds, oob error correct for user-specified function',
 code = {

  tree_seeds = sample.int(n = 50000, size = 100)
  bad_tree_seeds <- c(1,2,3)

  expect_error(
   orsf(data = pbc_orsf,
        formula = time+status~.-id,
        n_tree = 100,
        mtry = 2,
        tree_seeds = bad_tree_seeds),
   regexp = 'the number of trees'
  )

  fit_1 <- orsf(data = pbc_orsf,
                formula = time+status~.-id,
                n_tree = 100,
                mtry = 2,
                tree_seeds = tree_seeds)

  fit_2 <- orsf(data = pbc_orsf,
                formula = time+status~.-id,
                n_tree = 100,
                mtry = 6,
                tree_seeds = tree_seeds)

  expect_equal(fit_1$forest$rows_oobag,
               fit_2$forest$rows_oobag)

  fit_3 <- orsf(data = pbc_orsf,
                formula = time+status~.-id,
                n_tree = 100,
                mtry = 6,
                oobag_fun = oobag_c_survival,
                tree_seeds = tree_seeds)

  expect_equal(
   oobag_c_survival(
    y_mat = as.matrix(pbc_orsf[,c("time", "status")]),
    w_vec = rep(1, nrow(pbc_orsf)),
    s_vec = fit_3$pred_oobag
   ),
   as.numeric(fit_3$eval_oobag$stat_values)
  )

 }
)


if(Sys.getenv("run_all_aorsf_tests") == 'yes'){

 test_that(
  desc = 'orsf_time_to_train is reasonable at approximating time to train',
  code = {

   # testing the seed behavior when no_fit is TRUE. You should get the same
   # forest whether you train with orsf() or with orsf_train().

   for(.n_tree in c(100, 250, 1000)){

    object <- orsf(pbc_orsf, Surv(time, status) ~ . - id,
                   n_tree = .n_tree, no_fit = TRUE,
                   importance = 'anova')
    set.seed(89)
    time_estimated <- orsf_time_to_train(object, n_tree_subset = 50)

    set.seed(89)
    time_true_start <- Sys.time()
    fit_orsf_3 <- orsf_train(object)
    time_true_stop <- Sys.time()

    time_true <- time_true_stop - time_true_start

    diff_abs <- abs(as.numeric(time_true - time_estimated))
    diff_rel <- diff_abs / as.numeric(time_true)

    # expect the difference between estimated and true time is < 5 second.
    expect_lt(diff_abs, 5)
    # expect that the difference is not greater than 5x the
    # magnitude of the actual time it took to fit the forest
    expect_lt(diff_rel, 5)

   }
  }
 )

}


test_that(
 desc = 'orsf_train does not accept bad inputs',
 code = {

  expect_error(orsf_train(object = fit_orsf), regexp = 'been trained')

  fit_nodat <- orsf(pbc_orsf,
                    Surv(time, status) ~ . - id,
                    n_tree = 2,
                    no_fit = TRUE,
                    attach_data = FALSE)

  expect_error(orsf_train(object = fit_nodat),
               regexp = 'training data attached')

 }
)

test_that(
 desc = "results are similar after adding trivial noise",
 code = {

  expect_true(
   abs(fit_orsf$eval_oobag$stat_values - fit_orsf_noise$eval_oobag$stat_values) < 0.01
  )

  expect_true(
   mean(abs(fit_orsf$pred_oobag-fit_orsf_noise$pred_oobag)) < 0.1
  )

 }

)

# test_that(
#  desc = 'orsf_fit objects can be saved and loaded with saveRDS and readRDS',
#  code = {
#
#   fil <- tempfile("fit_orsf", fileext = ".rds")
#
#   ## save a single object to file
#   saveRDS(fit_orsf, fil)
#   ## restore it under a different name
#   fit_orsf_read_in <- readRDS(fil)
#
#   # NULL these attributes because they are functions
#   # the env of functions in fit_orsf_read_in will not be identical to the env
#   # of functions in fit_orsf. Everything else should be identical.
#
#   attr(fit_orsf, 'f_beta') <- NULL
#   attr(fit_orsf_read_in, 'f_beta') <- NULL
#
#   attr(fit_orsf, 'f_oobag_eval') <- NULL
#   attr(fit_orsf_read_in, 'f_oobag_eval') <- NULL
#
#   expect_equal(fit_orsf, fit_orsf_read_in)
#
#   p1=predict(fit_orsf,
#              new_data = fit_orsf$data,
#              pred_horizon = 1000)
#
#   p2=predict(fit_orsf_read_in,
#              new_data = fit_orsf_read_in$data,
#              pred_horizon = 1000)
#
#   expect_equal(p1, p2)
#
#  }
# )



#' @srrstats {ML7.9} *Explicitly compare all possible combinations in categorical differences in model architecture, such as different model architectures with same optimization algorithms, same model architectures with different optimization algorithms, and differences in both.*


test_that(
 desc = 'orsf() runs as intended for many valid inputs',
 code = {

  #' @srrstats {ML7.9a} *form combinations of inputs using `expand.grid()`.*
  inputs <- expand.grid(
   data_format = c('plain', 'tibble', 'data.table'),
   n_tree = 1,
   n_split = 1,
   n_retry = 0,
   mtry = 3,
   leaf_min_events = 5,
   leaf_min_obs = c(10),
   split_rule = c("logrank", "cstat"),
   split_min_events = 5,
   split_min_obs = 15,
   oobag_pred_type = c('none', 'risk', 'surv', 'chf', 'mort'),
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
                     'cph' = orsf_control_cph(),
                     'net' = orsf_control_net(),
                     'custom' = orsf_control_custom(beta_fun = f_pca))

   fit <- orsf(data = data_fun(pbc_orsf),
               formula = time + status ~ . - id,
               control = control,
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

   expect_s3_class(fit, class = 'orsf_fit')
   expect_equal(get_n_tree(fit), inputs$n_tree[i])
   expect_equal(get_n_split(fit), inputs$n_split[i])
   expect_equal(get_n_retry(fit), inputs$n_retry[i])
   expect_equal(get_mtry(fit), inputs$mtry[i])
   expect_equal(get_leaf_min_events(fit), inputs$leaf_min_events[i])
   expect_equal(get_leaf_min_obs(fit), inputs$leaf_min_obs[i])
   expect_equal(get_split_min_events(fit), inputs$split_min_events[i])
   expect_equal(get_split_min_obs(fit), inputs$split_min_obs[i])
   expect_equal(fit$pred_horizon, pred_horizon)

   expect_length(fit$forest$rows_oobag, n = get_n_tree(fit))

   if(inputs$oobag_pred_type[i] != 'none'){

    expect_length(fit$eval_oobag$stat_values, length(pred_horizon))
    expect_equal(nrow(fit$pred_oobag), get_n_obs(fit))

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
 desc = 'if oobag time is unspecified, pred horizon = median(time)',
 code = {

  fit_1 <- orsf(data = pbc_orsf,
                formula = time + status ~ . - id,
                n_tree = 1)

  fit_2 <- orsf(data = pbc_orsf,
                formula = time + status ~ . - id,
                n_tree = 1,
                oobag_pred_type = 'none')

  expect_equal(fit_1$pred_horizon, fit_2$pred_horizon)

 }
)

pbc_temp <- pbc_orsf
pbc_temp$list_col <- list(list(a=1))

#' @srrstats {G2.12} *pre-processing identifies list columns and throws informative error*

test_that(
 desc = 'list columns are not allowed',
 code = {
  expect_error(
   orsf(pbc_temp, time + status ~ . - id),
   regexp = '<list_col>'
  )
 }
)

set.seed(329)

fit_unwtd <- orsf(pbc_orsf, Surv(time, status) ~ . - id)

fit_wtd <- orsf(pbc_orsf, Surv(time, status) ~ . - id,
                weights = rep(2, nrow(pbc_orsf)))

test_that(
 desc = 'weights work as intended',
 code = {

  # using weights should make the trees much deeper:
  expect_gt(get_n_leaves_mean(fit_wtd),
            get_n_leaves_mean(fit_unwtd))

 }
)

test_that(
 desc = "lists can be plugged into orsf",
 code = {

  pbc_list <- as.list(pbc_orsf)
  pbc_list_bad <- pbc_list
  pbc_list_bad$trt <- pbc_list_bad$trt[1:3]
  pbc_list_bad$age <- pbc_list_bad$age[1:5]

  # only run locally - I don't want to list recipes in suggests
  # recipe <- recipes::recipe(pbc_orsf, formula = time + status ~ .) %>%
  #  recipes::step_rm(id) %>%
  #  recipes::step_scale(recipes::all_numeric_predictors())
  #
  # recipe_prepped <- recipes::prep(recipe)
  #
  # fit_recipe <- orsf(recipe_prepped, Surv(time, status) ~ .)
  #
  # expect_s3_class(fit_recipe, 'orsf_fit')

  fit_list <- orsf(pbc_list, Surv(time, status) ~ .)

  expect_s3_class(fit_list, 'orsf_fit')

  expect_error(
   orsf(pbc_list_bad, Surv(time, status) ~ .),
   regexp = 'unable to cast data'
  )

 }
)

# high pred horizon
# TODO: move this to test file for summarize
# test_that(
#  desc = 'higher pred horizon is not allowed for summary',
#  code = {
#
#   fit_bad_oob_horizon <- orsf(time + status ~ ., data = pbc_orsf,
#                               oobag_pred_horizon = 7000)
#
#   expect_error(orsf_summarize_uni(fit_bad_oob_horizon),
#                regexp = 'prediction horizon')
#
#  }
# )


# Similar to obliqueRSF?
# suppressPackageStartupMessages({
#  library(obliqueRSF)
# })
#
# set.seed(50)
#
# fit_aorsf <- orsf(pbc_orsf,
#                   formula = Surv(time, status) ~ . - id,
#                   n_tree = 100)
# fit_obliqueRSF <- ORSF(pbc_orsf, ntree = 100, verbose = FALSE)
#
#
# risk_aorsf <- predict(fit_aorsf, new_data = pbc_orsf, pred_horizon = 3500)
# risk_obliqueRSF <- 1-predict(fit_obliqueRSF, newdata = pbc_orsf, times = 3500)
#
# cor(risk_obliqueRSF, risk_aorsf)
# plot(risk_obliqueRSF, risk_aorsf)


