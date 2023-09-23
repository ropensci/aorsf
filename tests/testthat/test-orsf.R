
library(survival) # for Surv

# misc functions used for tests ----

no_miss_list <- function(l){

 sapply(l, function(x){

  if(is.list(x)) {return(no_miss_list(x))}

  any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x))

 })

}

add_noise <- function(x, eps = .Machine$double.eps){

 noise <- rnorm(length(x), mean = 0, sd = eps/2)
 noise <- pmin(noise, eps)
 noise <- pmax(noise, -eps)

 x + noise

}

change_scale <- function(x, mult_by = 1/2){
 x * mult_by
}

# begin tests -----

#' @srrstats {G5.0} *tests use the PBC data, a standard set that has been widely studied and disseminated in other R package (e.g., survival and randomForestSRC)*

# catch bad inputs, give informative error

pbc_temp <- pbc_orsf
pbc_temp$id <- factor(pbc_temp$id)
pbc_temp$status <- pbc_temp$status+1


f1 <- Surv(time, status) ~ unknown_variable + bili
# dropped test - see https://github.com/mlr-org/mlr3extralearners/issues/259
# f2 <- Surv(time, status) ~ bili
f3 <- Surv(time, status) ~ bili + factor(hepato)
f4 <- Surv(time, status) ~ bili * ascites
f5 <- Surv(time, status) ~ bili + id
f6 <- Surv(time, not_right) ~ .
f7 <- Surv(not_right, status) ~ .
f8 <- Surv(start, time, status) ~ .
f9 <- Surv(status, time) ~ . - id
f10 <- Surv(time, time) ~ . - id
f11 <- Surv(time, id) ~ . -id
f12 <- Surv(time, status) ~ . -id
f13 <- ~ .
f14 <- status + time ~ . - id
f15 <- time + status ~ id + bili

#' @srrstats {G5.2} *Appropriate error behaviour is explicitly demonstrated through tests.*
#' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*
test_that(
 desc = 'formula inputs are vetted',
 code = {

  expect_error(orsf(pbc_temp, f1), 'not found in data')
  # # dropped - see https://github.com/mlr-org/mlr3extralearners/issues/259
  # expect_warning(orsf(pbc_temp, f2), 'at least 2 predictors')
  expect_error(orsf(pbc_temp, f3), 'unrecognized')
  expect_error(orsf(pbc_temp, f4), 'unrecognized')
  expect_error(orsf(pbc_temp, f5), 'id variable?')
  expect_error(orsf(pbc_temp, f6), 'not_right')
  expect_error(orsf(pbc_temp, f7), 'not_right')
  expect_error(orsf(pbc_temp, f8), 'must have two variables')
  expect_error(orsf(pbc_temp, f9), 'Did you enter')
  expect_error(orsf(pbc_temp, f10), 'must have two variables')
  expect_error(orsf(pbc_temp, f11), 'detected >1 event type')
  expect_error(orsf(pbc_temp, f13), 'must be two sided')
  expect_error(orsf(pbc_temp, f14), 'Did you enter')
  expect_error(orsf(pbc_temp, f15), "as many levels as there are rows")

 }
)

test_that(
 desc = 'long formulas with repetition are allowed',
 code = {

  x_vars <- c(
   "trt",
   "age",
   "sex",
   "ascites",
   "hepato",
   "spiders",
   "edema",
   "bili",
   "chol",
   "albumin",
   "copper",
   "alk.phos",
   "ast",
   "trig",
   "platelet",
   "protime",
   "stage"
  )

  long_rhs <- paste(x_vars, collapse = ' + ')

  long_rhs <- rep(long_rhs, 15)

  long_rhs <- paste(long_rhs, collapse = ' + ')

  f_long <- as.formula(paste("time + status ~", long_rhs))

  fit_long <- orsf(formula = f_long, pbc_temp, n_tree = 10)

  # fits the orsf as expected
  expect_s3_class(fit_long, 'orsf_fit')
  # keeps unique names
  expect_equal(x_vars, get_names_x(fit_long))

 }
)

# should get the same forest, whether status is 1/2 or 0/1 or a surv object

pbc_surv <- Surv(pbc_temp$time, pbc_temp$status)
pbc_surv_data <- cbind(pbc_temp, surv_object=pbc_surv)

fit_surv <- orsf(pbc_surv_data,
                 formula = surv_object ~ . - id - time - status,
                 n_tree = 10,
                 tree_seed = 1:10)

fit_surv_untrained <- orsf(pbc_surv_data,
                           formula = surv_object ~ . - id - time - status,
                           n_tree = 10,
                           tree_seed = 1:10,
                           no_fit = TRUE)

fit_surv_trained <- orsf_train(fit_surv_untrained)

fit_12 <- orsf(pbc_temp,
               formula = Surv(time, status) ~ . -id,
               n_tree = 10,
               tree_seeds = 1:10)

fit_01 <- orsf(pbc_orsf,
               formula = time + status ~ . -id,
               n_tree = 10,
               tree_seeds = 1:10)


test_that(
 desc = 'New status, same forest',
 code = {
  expect_identical(fit_12$forest, fit_01$forest)
  expect_identical(fit_surv$forest, fit_01$forest)
  expect_identical(fit_surv_trained$forest, fit_01$forest)
 }
)


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

  pbc_temp$date_var <- Sys.Date()
  expect_error(orsf(pbc_temp, f), 'unsupported type')

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
      n_tree = 100,
      tree_seeds = 1:100)

fit_orsf_2 <-
 orsf(pbc_orsf, Surv(time, status) ~ . - id,
      n_thread = 5,
      n_tree = 100,
      tree_seeds = 1:100)

fit_orsf_noise <-
 orsf(pbc_noise, Surv(time, status) ~ . - id,
      n_tree = 100,
      tree_seeds = 1:100)

fit_orsf_scale <-
 orsf(pbc_scale, Surv(time, status) ~ . - id,
      n_tree = 100,
      tree_seeds = 1:100)

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
   mean(abs(fit_orsf$pred_oobag - fit_orsf_scale$pred_oobag)),
   0.1
  )

  expect_lt(
   mean(abs(fit_orsf$pred_oobag - fit_orsf_2$pred_oobag)),
   0.1
  )

  expect_lt(
   mean(abs(fit_orsf$pred_oobag - fit_orsf_noise$pred_oobag)),
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
                 n_tree = 100,
                 tree_seeds = 1:100,
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
 desc = 'orsf() runs as intended across numerous possible architectures',
 code = {

  #' @srrstats {ML7.9a} *form combinations of inputs using `expand.grid()`.*
  inputs <- expand.grid(
   data_format = c('plain', 'tibble', 'data.table'),
   n_tree = 1,
   n_split = 1,
   n_retry = 0,
   mtry = 3,
   leaf_min_events = 1,
   leaf_min_obs = c(5, 10),
   split_min_events = 5,
   split_min_obs = 15,
   oobag_pred_type = c('none', 'risk', 'surv', 'chf'),
   oobag_pred_horizon = c(1000),
   stringsAsFactors = FALSE
  )

  for(i in seq(nrow(inputs))){

   data_fun <- switch(
    as.character(inputs$data_format[i]),
    'plain' = function(x) x,
    'tibble' = tibble::as_tibble,
    'data.table' = as.data.table
   )

   fit_cph <- orsf(data = data_fun(pbc_orsf),
                   formula = time + status ~ . - id,
                   control = orsf_control_cph(),
                   n_tree = inputs$n_tree[i],
                   n_split = inputs$n_split[i],
                   n_retry = inputs$n_retry[i],
                   mtry = inputs$mtry[i],
                   leaf_min_events = inputs$leaf_min_events[i],
                   leaf_min_obs = inputs$leaf_min_obs[i],
                   split_min_events = inputs$split_min_events[i],
                   split_min_obs = inputs$split_min_obs[i],
                   oobag_pred_type = inputs$oobag_pred_type[i],
                   oobag_pred_horizon = inputs$oobag_pred_horizon[i])

   expect_s3_class(fit_cph, class = 'orsf_fit')
   expect_equal(get_n_tree(fit_cph), inputs$n_tree[i])
   expect_equal(get_n_split(fit_cph), inputs$n_split[i])
   expect_equal(get_n_retry(fit_cph), inputs$n_retry[i])
   expect_equal(get_mtry(fit_cph), inputs$mtry[i])
   expect_equal(get_leaf_min_events(fit_cph), inputs$leaf_min_events[i])
   expect_equal(get_leaf_min_obs(fit_cph), inputs$leaf_min_obs[i])
   expect_equal(get_split_min_events(fit_cph), inputs$split_min_events[i])
   expect_equal(get_split_min_obs(fit_cph), inputs$split_min_obs[i])
   expect_equal(fit_cph$pred_horizon, inputs$oobag_pred_horizon[i])

   expect_length(fit_cph$forest$rows_oobag, n = get_n_tree(fit_cph))

   if(inputs$oobag_pred_type[i] != 'none'){
    expect_length(fit_cph$eval_oobag$stat_values, 1)
    expect_equal(nrow(fit_cph$pred_oobag), get_n_obs(fit_cph))
   } else {
    expect_equal(dim(fit_cph$eval_oobag$stat_values), c(0, 0))
   }

   fit_net <- orsf(data = pbc_orsf,
                   formula = time + status ~ . - id,
                   control = orsf_control_net(),
                   n_tree = 1,
                   n_split = inputs$n_split[i],
                   n_retry = inputs$n_retry[i],
                   mtry = inputs$mtry[i],
                   leaf_min_events = inputs$leaf_min_events[i],
                   leaf_min_obs = inputs$leaf_min_obs[i],
                   split_min_events = inputs$split_min_events[i],
                   split_min_obs = inputs$split_min_obs[i],
                   oobag_pred_type = inputs$oobag_pred_type[i],
                   oobag_pred_horizon = inputs$oobag_pred_horizon[i])


   expect_s3_class(fit_net, class = 'orsf_fit')
   expect_equal(get_n_tree(fit_net), inputs$n_tree[i])
   expect_equal(get_n_split(fit_net), inputs$n_split[i])
   expect_equal(get_n_retry(fit_net), inputs$n_retry[i])
   expect_equal(get_mtry(fit_net), inputs$mtry[i])
   expect_equal(get_leaf_min_events(fit_net), inputs$leaf_min_events[i])
   expect_equal(get_leaf_min_obs(fit_net), inputs$leaf_min_obs[i])
   expect_equal(get_split_min_events(fit_net), inputs$split_min_events[i])
   expect_equal(get_split_min_obs(fit_net), inputs$split_min_obs[i])
   expect_equal(fit_net$pred_horizon, inputs$oobag_pred_horizon[i])

   expect_length(fit_cph$forest$rows_oobag, n = get_n_tree(fit_cph))

   if(inputs$oobag_pred_type[i] != 'none'){
    expect_length(fit_net$eval_oobag$stat_values, 1)
    expect_equal(nrow(fit_net$pred_oobag), get_n_obs(fit_net))
   } else {
    expect_equal(dim(fit_net$eval_oobag$stat_values), c(0, 0))
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


