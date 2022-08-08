oobag_c_harrell <- function(y_mat, s_vec){

 sorted <- order(y_mat[, 1], -y_mat[, 2])

 y_mat <- y_mat[sorted, ]
 s_vec <- s_vec[sorted]

 time = y_mat[, 1]
 status = y_mat[, 2]
 events = which(status == 1)

 k = nrow(y_mat)

 total <- 0
 concordant <- 0

 for(i in events){

  if(i+1 <= k){

   for(j in seq(i+1, k)){

    if(time[j] > time[i]){

     total <- total + 1

     if(s_vec[j] > s_vec[i]){

      concordant <- concordant + 1

     } else if (s_vec[j] == s_vec[i]){

      concordant <- concordant + 0.5

     }

    }

   }

  }

 }

 concordant / total

}

#' @srrstats {G5.0} *tests use the PBC data, a standard set that has been widely studied and disseminated in other R package (e.g., survival and randomForestSRC)*

# catch bad inputs, give informative error

pbc_temp <- pbc_orsf
pbc_temp$id <- factor(pbc_temp$id)
pbc_temp$status <- pbc_temp$status+1


f1 <- Surv(time, status) ~ unknown_variable + bili
f2 <- Surv(time, status) ~ id
f3 <- Surv(time, status) ~ bili + factor(hepato)
f4 <- Surv(time, status) ~ bili * ascites
f5 <- Surv(time, status) ~ bili + id
f6 <- Surv(time, not_right) ~ .
f7 <- Surv(not_right, status) ~ .
f8 <- Surv(start, time, status) ~ .
f9 <- Surv(status, time) ~ . - id
f10 <- Surv(time, time) ~ . - id
f11 <- Surv(time, hepato) ~ . -id
f12 <- Surv(time, status) ~ . -id
f13 <- ~ .
f14 <- status + time ~ . - id

#' @srrstats {G5.2} *Appropriate error behaviour is explicitly demonstrated through tests.*
#' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*
test_that(
 desc = 'formula inputs are vetted',
 code = {

  expect_error(orsf(pbc_temp, f1), 'not found in data')
  expect_error(orsf(pbc_temp, f2), 'at least 2 predictors')
  expect_error(orsf(pbc_temp, f3), 'unrecognized')
  expect_error(orsf(pbc_temp, f4), 'unrecognized')
  expect_error(orsf(pbc_temp, f5), 'id variable?')
  expect_error(orsf(pbc_temp, f6), 'not_right')
  expect_error(orsf(pbc_temp, f7), 'not_right')
  expect_error(orsf(pbc_temp, f8), 'must have two variables')
  expect_error(orsf(pbc_temp, f9), 'should contain values of 0 and 1')
  expect_error(orsf(pbc_temp, f10), 'must have two variables')
  expect_error(orsf(pbc_temp, f11), 'should have type')
  expect_error(orsf(pbc_temp, f12), 'should contain values of 0 and 1')
  expect_error(orsf(pbc_temp, f13), 'must be two sided')
  expect_error(orsf(pbc_temp, f14), 'should contain values of 0 and 1')

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

  pbc_temp$date_var <- Sys.Date()
  expect_error(orsf(pbc_temp, f), 'unsupported type')

 }
)

test_that(
 desc = 'target_df too high is caught',
 code = {

  cntrl <- orsf_control_net(df_target = 10)
  expect_error(orsf(pbc_orsf, f, cntrl), 'must be <= mtry')

 }
)

test_that(
 desc = 'orsf runs with data.table and with net control',
 code = {

  expect_s3_class(orsf(as.data.table(pbc_orsf), f, n_tree = 1), 'aorsf')

  expect_s3_class(orsf(as.data.table(pbc_orsf), f,
                       control = orsf_control_net(),
                       n_tree = 1), 'aorsf')
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

fit_with_vi <- orsf(data = pbc_orsf,
                    formula = Surv(time, status) ~ . - id,
                    importance = 'negate',
                    n_tree = 50)

fit_no_vi <- orsf(data = pbc_orsf,
                  formula = Surv(time, status) ~ . - id,
                  importance = 'none',
                  n_tree = 50)

no_miss_list <- function(l){

 sapply(l, function(x){

  if(is.list(x)) {return(no_miss_list(x))}

  any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x))

 })

}

#' @srrstats {G5.3} *Explicit test expected to return objects containing no missing (`NA`) or undefined (`NaN`, `Inf`) values are explicitly tested.*

test_that(
 "output contains no missing values",
 code = {

  miss_check_no_vi <- sapply(fit_no_vi, no_miss_list)
  miss_check_with_vi <- sapply(fit_with_vi, no_miss_list)

  for(i in seq_along(miss_check_no_vi)){
   if(!is_empty(miss_check_no_vi[[i]]))
    expect_true(sum(miss_check_no_vi[[i]]) == 0)
  }

  for(i in seq_along(miss_check_with_vi)){
   if(!is_empty(miss_check_with_vi[[i]]))
    expect_true(sum(miss_check_with_vi[[i]]) == 0)
  }


 }
)

cstat_bcj <- function(y_mat, s_vec){

 sorted <- order( y_mat[, 1], -y_mat[, 2])
 oobag_c_harrell_testthat(y_mat[sorted, ], s_vec[sorted, ])

}

test_that(
 desc = 'oobag error is reproducible from an aorsf object',
 code = {

  y_mat <- as.matrix(fit_no_vi$data[, c('time', 'status')])
  s_vec <- fit_no_vi$surv_oobag

  tt <- survival::concordancefit(
   y = survival::Surv(pbc_orsf$time, pbc_orsf$status),
   x = fit_no_vi$surv_oobag
  )

  denom <- sum(tt$count[c('concordant',
                          'discordant',
                          'tied.y')])

  target <- as.numeric(tt$concordance)

  bcj <- cstat_bcj(y_mat, s_vec)

  expect_equal(
   bcj,
   as.numeric(fit_no_vi$eval_oobag$stat_values)
  )

  # cstat_bcj close enough to cstat from survival
  expect_lt(abs(target - bcj), 0.001)

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
        control = orsf_control_cph(iter_max = 50, eps = 1),
        Surv(time, status) ~ . -id,
        n_tree = 150)
  )

  time_large <- system.time(
   orsf(pbc_orsf,
        control = orsf_control_cph(iter_max = 50, eps = 1e-10),
        Surv(time, status) ~ . -id,
        n_tree = 150)
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
 desc = "data with missing values are rejected",
 code = {
  expect_error(orsf(pbc_temp, time + status ~ . - id),
               'missing values')
 }
)

#' @srrstats {G5.9} **Noise susceptibility tests**
#' @srrstats {G5.9a} *Adding trivial noise to data does not meaningfully change results*
#' @srrstats {G5.9b} *Running under different random seeds gives identifal results*

add_noise <- function(x, eps = .Machine$double.eps){
 x + rnorm(length(x), mean = 0, sd = eps)
}

change_scale <- function(x, mult_by = 10){
 x * mult_by
}

pbc_noise <- pbc_orsf
pbc_scale <- pbc_orsf

vars <- c('bili', 'chol', 'albumin', 'copper', 'alk.phos', 'ast')

pbc_noise[, vars] <- sapply(pbc_noise[, vars], add_noise)

pbc_scale[, vars] <- sapply(pbc_scale[, vars], change_scale)

set.seed(329)

fit_orsf <- orsf(pbc_orsf,
                 Surv(time, status) ~ . - id,
                 n_tree = 10)

set.seed(329)

fit_orsf_2 <- orsf(pbc_orsf,
                   Surv(time, status) ~ . - id,
                   n_tree = 10)

set.seed(329)

fit_orsf_noise <- orsf(pbc_noise,
                       Surv(time, status) ~ . - id,
                       n_tree = 10)

set.seed(329)

fit_orsf_scale <- orsf(pbc_scale,
                       Surv(time, status) ~ . - id,
                       n_tree = 10)

#' @srrstats {ML7.1} *Demonstrate effect of numeric scaling of input data.*
test_that(
 desc = 'scaling/noising inputs does not impact model behavior',
 code = {

  expect_equal(fit_orsf$eval_oobag$stat_values,
               fit_orsf_scale$eval_oobag$stat_values)

  expect_equal(fit_orsf$eval_oobag$stat_values,
               fit_orsf_2$eval_oobag$stat_values)

  expect_equal(fit_orsf$eval_oobag$stat_values,
               fit_orsf_noise$eval_oobag$stat_values)

  # expect_equal(fit_orsf$surv_oobag,
  #              fit_orsf_scale$surv_oobag)
  #
  #   expect_equal(fit_orsf$surv_oobag,
  #                fit_orsf_2$surv_oobag)
  #
  #   expect_equal(fit_orsf$surv_oobag,
  #                fit_orsf_noise$surv_oobag)

  expect_equal(fit_orsf$forest[[1]]$leaf_nodes,
               fit_orsf_2$forest[[1]]$leaf_nodes)

  expect_equal(fit_orsf$forest[[1]]$leaf_nodes,
               fit_orsf_scale$forest[[1]]$leaf_nodes)

  expect_equal(fit_orsf$forest[[1]]$leaf_nodes,
               fit_orsf_noise$forest[[1]]$leaf_nodes)

 }
)

# testing the seed behavior when no_fit is TRUE. You should get the same
# forest whether you train with orsf() or with orsf_train().


object <- orsf(pbc_orsf, Surv(time, status) ~ . - id,
               n_tree = 10, no_fit = TRUE)

set.seed(329)

fit_orsf_3 <- orsf_train(object)

test_that(
 desc = 'results are identical if a forest is fitted under the same random seed',
 code = {

  # testing a subset of trees for identical betas

  for(i in seq(get_n_tree(fit_orsf))){
   expect_equal(
    object = fit_orsf$forest[[i]]$betas,
    expected = fit_orsf_2$forest[[i]]$betas
   )
   expect_equal(
    object = fit_orsf$forest[[i]]$betas,
    expected = fit_orsf_3$forest[[i]]$betas
   )
  }

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

  tree_seeds = sample.int(n = 50000, size = 10)
  bad_tree_seeds <- c(1,2,3)

  expect_error(
   orsf(data = pbc_orsf,
        formula = time+status~.-id,
        n_tree = 10,
        mtry = 2,
        tree_seeds = bad_tree_seeds),
   regexp = 'the number of trees'
  )

  fit_1 <- orsf(data = pbc_orsf,
                formula = time+status~.-id,
                n_tree = 10,
                mtry = 2,
                tree_seeds = tree_seeds)

  fit_2 <- orsf(data = pbc_orsf,
                formula = time+status~.-id,
                n_tree = 10,
                mtry = 6,
                tree_seeds = tree_seeds)

  fit_3 <- orsf(data = pbc_orsf,
                formula = time+status~.-id,
                n_tree = 10,
                mtry = 6,
                oobag_fun = oobag_c_harrell,
                tree_seeds = tree_seeds)

  expect_equal(
   fit_2$eval_oobag$stat_values,
   fit_3$eval_oobag$stat_values
  )

  for(i in seq(get_n_tree(fit_2))){

   expect_equal(fit_1$forest[[i]]$rows_oobag,
                fit_2$forest[[i]]$rows_oobag)

  }
 }
)



test_that(
 desc = 'orsf_time_to_train is reasonable at approximating time to train',
 code = {

  # testing the seed behavior when no_fit is TRUE. You should get the same
  # forest whether you train with orsf() or with orsf_train().

  for(.n_tree in c(100, 250, 1000, 2500)){

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
   max(abs(fit_orsf$surv_oobag-fit_orsf_noise$surv_oobag)) < 0.1
  )

 }

)

test_that(
 desc = 'aorsf objects can be saved and loaded with saveRDS and readRDS',
 code = {

  fil <- tempfile("fit_orsf", fileext = ".rds")

  ## save a single object to file
  saveRDS(fit_orsf, fil)
  ## restore it under a different name
  fit_orsf_read_in <- readRDS(fil)

  # NULL these attributes because they are functions
  # the env of functions in fit_orsf_read_in will not be identical to the env
  # of functions in fit_orsf. Everything else should be identical.

  attr(fit_orsf, 'f_beta') <- NULL
  attr(fit_orsf_read_in, 'f_beta') <- NULL

  attr(fit_orsf, 'f_oobag_eval') <- NULL
  attr(fit_orsf_read_in, 'f_oobag_eval') <- NULL

  expect_equal(fit_orsf, fit_orsf_read_in)

  p1=predict(fit_orsf,
             new_data = fit_orsf$data,
             pred_horizon = 1000)

  p2=predict(fit_orsf_read_in,
             new_data = fit_orsf_read_in$data,
             pred_horizon = 1000)

  expect_equal(p1, p2)

 }
)



#' @srrstats {ML7.9} *Explicitly compare all possible combinations in categorical differences in model architecture, such as different model architectures with same optimization algorithms, same model architectures with different optimization algorithms, and differences in both.*


test_that(
 desc = 'orsf() runs as intended across numerous possible architectures',
 code = {

  #' @srrstats {ML7.9a} *form combinations of inputs using `expand.grid()`.*
  inputs <- expand.grid(
   n_tree = 1,
   n_split = 1,
   n_retry = c(0, 3),
   mtry = c(2, 12),
   leaf_min_events = c(1, 3),
   leaf_min_obs = c(5, 10),
   split_min_events = c(6, 9),
   split_min_obs = 15,
   oobag_pred = c(TRUE, FALSE),
   oobag_pred_horizon = c(2000, 4000)
  )

  for(i in seq(nrow(inputs))){

   fit_cph <- orsf(data = pbc_orsf,
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
                   oobag_pred = inputs$oobag_pred[i],
                   oobag_pred_horizon = inputs$oobag_pred_horizon[i])

   expect_s3_class(fit_cph, class = 'aorsf')
   expect_equal(get_n_tree(fit_cph), inputs$n_tree[i])
   expect_equal(get_n_split(fit_cph), inputs$n_split[i])
   expect_equal(get_n_retry(fit_cph), inputs$n_retry[i])
   expect_equal(get_mtry(fit_cph), inputs$mtry[i])
   expect_equal(get_leaf_min_events(fit_cph), inputs$leaf_min_events[i])
   expect_equal(get_leaf_min_obs(fit_cph), inputs$leaf_min_obs[i])
   expect_equal(get_split_min_events(fit_cph), inputs$split_min_events[i])
   expect_equal(get_split_min_obs(fit_cph), inputs$split_min_obs[i])
   expect_equal(get_oobag_pred(fit_cph), inputs$oobag_pred[i])
   expect_equal(fit_cph$pred_horizon, inputs$oobag_pred_horizon[i])

   expect_length(fit_cph$forest, n = get_n_tree(fit_cph))

   if(inputs$oobag_pred[i]){
    expect_length(fit_cph$eval_oobag$stat_values, 1)
    expect_equal(nrow(fit_cph$surv_oobag), get_n_obs(fit_cph))
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
                   oobag_pred = inputs$oobag_pred[i],
                   oobag_pred_horizon = inputs$oobag_pred_horizon[i])

   expect_s3_class(fit_net, class = 'aorsf')
   expect_equal(get_n_tree(fit_net), inputs$n_tree[i])
   expect_equal(get_n_split(fit_net), inputs$n_split[i])
   expect_equal(get_n_retry(fit_net), inputs$n_retry[i])
   expect_equal(get_mtry(fit_net), inputs$mtry[i])
   expect_equal(get_leaf_min_events(fit_net), inputs$leaf_min_events[i])
   expect_equal(get_leaf_min_obs(fit_net), inputs$leaf_min_obs[i])
   expect_equal(get_split_min_events(fit_net), inputs$split_min_events[i])
   expect_equal(get_split_min_obs(fit_net), inputs$split_min_obs[i])
   expect_equal(get_oobag_pred(fit_net), inputs$oobag_pred[i])
   expect_equal(fit_net$pred_horizon, inputs$oobag_pred_horizon[i])

   if(inputs$oobag_pred[i]){
    expect_length(fit_net$eval_oobag$stat_values, 1)
    expect_equal(nrow(fit_net$surv_oobag), get_n_obs(fit_net))
   }

  }

  test_that(
   desc = 'if oobag time is unspecified, pred horizon = median(time)',
   code = {

    fit_1 <- orsf(data = pbc_orsf,
                  formula = time + status ~ . - id,
                  n_tree = 1,
                  oobag_pred = TRUE)

    fit_2 <- orsf(data = pbc_orsf,
                  formula = time + status ~ . - id,
                  n_tree = 1,
                  oobag_pred = FALSE)

    expect_equal(fit_1$pred_horizon, fit_2$pred_horizon)

   }
  )

 }

)

pbc_temp <- pbc_orsf
pbc_temp$list_col <- list(a=1)

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
                weights = pbc_orsf$id)

# using weights should make the trees much deeper:
expect_gt(get_n_leaves_mean(fit_wtd),
          get_n_leaves_mean(fit_unwtd))

# and in this case less accurate b/c the weights were random and extreme
expect_lt(
 fit_wtd$eval_oobag$stat_values,
 fit_unwtd$eval_oobag$stat_values
)




