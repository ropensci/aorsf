
test_preds_surv <- function(pred_type){

 n_train <- nrow(pbc_train)
 n_test <- nrow(pbc_test)

 pred_ncols_expect_agg <- switch(
  pred_type,
  risk = length(pred_horizon),
  surv = length(pred_horizon),
  chf  = length(pred_horizon),
  mort = 1,
  leaf = n_tree_test
 )

 dim_expect_agg <- list(
  oob = c(n_train, pred_ncols_expect_agg),
  new = c(n_test, pred_ncols_expect_agg)
 )

 dim_expect_raw <- c(n_test, n_tree_test)

 if(pred_type %in% c("chf", "risk", "surv"))
  dim_expect_raw <- c(dim_expect_raw, length(pred_horizon))

 fit <- orsf(formula = time + status  ~ . - id,
             data = pbc_train,
             oobag_pred_type = pred_type,
             n_tree = n_tree_test,
             oobag_pred_horizon = pred_horizon,
             tree_seeds = seeds_standard)

 if(pred_type %in% c("mort", "leaf")) pred_horizon <- NULL

 prd_agg <- predict(fit,
                    new_data = pbc_test,
                    pred_type = pred_type,
                    pred_horizon = pred_horizon,
                    n_thread = 1)

 prd_raw <- predict(fit,
                    new_data = pbc_test,
                    pred_aggregate = FALSE,
                    pred_type = pred_type,
                    pred_horizon = pred_horizon,
                    n_thread = 1)

 test_that(
  'No missing, nan, or infinite values in prediction output',
  code = {
   expect_false(any(is.na(prd_agg)))
   expect_false(any(is.nan(prd_agg)))
   expect_false(any(is.infinite(prd_agg)))
   expect_false(any(is.na(prd_raw)))
   expect_false(any(is.nan(prd_raw)))
   expect_false(any(is.infinite(prd_raw)))
  }
 )


 if(pred_type %in% c("risk", "surv")){

  test_that(
   desc = paste("predictions of type", pred_type, "are bounded"),
   code = {
    expect_true(all(prd_raw <= 1))
    expect_true(all(prd_raw >= 0))
   }
  )

 }

 if(pred_type == 'mort'){

  test_that(
   desc = "predictions are accurate",
   code = {

    surv_concord <- survival::concordance(
     survival::Surv(time, status) ~ prd_agg,
     data = pbc_test
    )

    mort_cstat <- 1 - surv_concord$concordance

    expect_true(mort_cstat > 0.60)

   }
  )

 }

 test_that(
  desc = paste(pred_type, "prediction dimensions match expectations"),
  code = {
   expect_equal(dim_expect_agg$oob, dim(fit$pred_oobag))
   expect_equal(dim_expect_agg$new, dim(prd_agg))
   expect_equal(dim_expect_raw, dim(prd_raw))
  }
 )

 test_that(
  desc = paste('thread stability for predictions of type', pred_type),
  code = {

   expect_equal(
    prd_agg,
    predict(fit,
            new_data = pbc_test,
            pred_type = pred_type,
            pred_horizon = pred_horizon,
            n_thread = 3)
   )

   expect_equal(
    prd_raw,
    predict(fit,
            new_data = pbc_test,
            pred_aggregate = FALSE,
            pred_type = pred_type,
            pred_horizon = pred_horizon,
            n_thread = 3)
   )

  }
 )

 list(fit = fit, prd_agg = prd_agg, prd_raw = prd_raw)

}

pred_horizon <- c(1000, 2500)

pred_objects_surv <- lapply(pred_types_surv, test_preds_surv)

test_that(
 desc = "prediction at time 0 is correct",
 code = {

  for(i in c("surv", "chf", "risk")){

   pred_t0 <- predict(fit_standard_pbc$fast,
                      new_data = pbc_test[1, ],
                      pred_type = i,
                      pred_horizon = 0)

   if(i %in% c("risk", "chf")) expect_equal(pred_t0, matrix(0))
   if(i %in% c("surv")) expect_equal(pred_t0, matrix(1))

  }
 }
)

test_that(
 desc = "leaf predictions aggregate same as raw",
 code = {
  expect_equal(pred_objects_surv$leaf$prd_raw,
               pred_objects_surv$leaf$prd_agg)
 }
)

test_that(
 desc = "unaggregated predictions can reproduce aggregated ones",
 code = {

  for(i in c("surv", "risk", "chf")){
   for(j in seq_along(pred_horizon)){
    expect_equal(
     pred_objects_surv[[i]]$prd_agg[, j],
     apply(pred_objects_surv[[i]]$prd_raw[, , j], 1, mean),
     tolerance = 1e-9
    )
   }
  }

  expect_equal(
   pred_objects_surv$mort$prd_agg,
   matrix(apply(pred_objects_surv$mort$prd_raw, 1, mean), ncol = 1)
  )

 }
)

test_that(
 desc = "same predictions from the forest regardless of oob type",
 code = {

  risk_preds <- lapply(
   pred_objects_surv,
   function(object){
    predict(object$fit,
            new_data = pbc_test,
            pred_horizon = 3500,
            pred_type = 'risk')
   }
  )

  for( i in seq(2, length(risk_preds))){
   expect_equal(risk_preds[[i]], risk_preds[[1]])
  }

 }
)

test_that(
 desc = 'predict is type stable',
 code = {
  for(i in seq_along(pred_objects_surv)){
   expect_true(is.array(pred_objects_surv[[i]]$prd_raw))
   expect_true(is.matrix(pred_objects_surv[[i]]$prd_agg))
  }
 }
)

# from here out we just test general predict() mechanics

fit <- fit_standard_pbc$fast

test_that(
 desc = "warnings served if pred_horizon is not needed",
 code = {

  expect_warning(
   predict(fit,
           new_data = pbc_orsf[1, ],
           pred_horizon = c(50, 500),
           pred_type = 'leaf'),
   regexp = 'does not impact predictions'
  )

  expect_warning(
   predict(fit,
           new_data = pbc_orsf[1, ],
           pred_horizon = c(50, 500),
           pred_type = 'mort'),
   regexp = 'does not impact predictions'
  )


 }
)

new_data <- pbc_test[1:10, ]

test_that(
 desc = 'pred_horizon automatically set to object$pred_horizon if needed',
 code = {
  expect_equal(
   predict(fit, new_data = new_data, pred_horizon = fit$pred_horizon),
   predict(fit, new_data = new_data)
  )
 }
)

test_that(
 desc = 'identical na_action = pass/fail/impute/omit if no missing data',
 code = {
  expect_equal(
   predict(fit, new_data = new_data, na_action = 'fail'),
   predict(fit, new_data = new_data, na_action = 'pass')
  )
  expect_equal(
   predict(fit, new_data = new_data, na_action = 'fail'),
   predict(fit, new_data = new_data, na_action = 'impute_meanmode')
  )
  expect_equal(
   predict(fit, new_data = new_data, na_action = 'fail'),
   predict(fit, new_data = new_data, na_action = 'omit')
  )
 }
)


test_that(
 desc = 'predictions computed for tibbles, and data.tables',
 code = {

  new_data_dt <- as.data.table(new_data)
  new_data_tbl <- tibble::as_tibble(new_data)

  for(pred_type in c("risk", "chf", "surv")){

   p1 <- predict(fit,
                 new_data = new_data,
                 pred_type = pred_type,
                 pred_horizon = c(1000, 2500))

   p1_dt <- predict(fit,
                    new_data = new_data_dt,
                    pred_type = pred_type,
                    pred_horizon = c(1000, 2500))

   p1_tbl <- predict(fit,
                     new_data = new_data_tbl,
                     pred_type = pred_type,
                     pred_horizon = c(1000, 2500))

   expect_equal(p1, p1_dt)
   expect_equal(p1, p1_tbl)

  }

  for(pred_type in c("mort", "leaf")){

   p1 <- predict(fit,
                 new_data = new_data,
                 pred_type = pred_type)

   p1_dt <- predict(fit,
                    new_data = new_data_dt,
                    pred_type = pred_type)

   p1_tbl <- predict(fit,
                     new_data = new_data_tbl,
                     pred_type = pred_type)

   expect_equal(p1, p1_dt)
   expect_equal(p1, p1_tbl)

  }

 }

)



test_that(
 desc = 'multi-time pred values independent of previous time',
 code = {

  for(pred_type in c("surv", "risk", "chf")){
   expect_equal(
    predict(fit,
            new_data = new_data,
            pred_type = pred_type,
            pred_horizon = c(500, 1500, 2000))[, 3],
    predict(fit,
            new_data = new_data,
            pred_type = pred_type,
            pred_horizon = c(1000, 2000))[, 2]
   )
  }

 }
)

test_that(
 desc = 'risk is inverse of survival',
 code = {
  p_risk <- predict(fit, new_data = new_data, pred_type = 'risk')
  p_surv <- predict(fit, new_data = new_data, pred_type = 'surv')
  expect_equal(p_risk, 1-p_surv, tolerance = 1e-9)
 }
)

test_that(
 desc = 'leaf predictions do not depend on other observations in the data',
 code = {

  for(pred_type in pred_types_surv){

   p_all <- predict(fit, new_data = new_data, pred_type = pred_type)

   for(i in seq(nrow(new_data))){
    p_1row <- predict(fit, new_data = new_data[i,], pred_type = pred_type)
    expect_equal(p_1row, p_all[i, , drop=FALSE])
   }

  }

 }
)

test_that(
 'leaf predictions do not depend on order of the data',
 code = {

  for(pred_type in pred_types_surv){

   p_before <- predict(fit,
                       new_data = new_data,
                       pred_type = pred_type)

   new_order <- sample(nrow(new_data), replace = F)

   p_after <- predict(fit,
                      new_data = new_data[new_order, ],
                      pred_type = pred_type)

   expect_equal(p_before[new_order, , drop = FALSE], p_after)

  }
 }
)

test_that(
 "mistakenly named inputs are caught",
 code = {

  expect_error(
   predict(fit, newdata = new_data, pred_horizon = 1000),
   regexp = 'newdata'
  )

  expect_error(
   predict(fit, newdata = new_data, horizon = 1000),
   regexp = 'horizon'
  )

  expect_error(
   predict(fit, newdata = new_data, horizon = 1000, type = 'risk'),
   regexp = 'type'
  )

  expect_error(
   predict(fit, OK = 'risk'),
   regexp = 'OK'
  )

 }
)

test_that(
 desc = 'Boundary case: empty new data throw an error',
 code = {

  expect_error(
   predict(fit, new_data = new_data[c(), ], pred_horizon = 1000),
   regexp = 'new data are empty'
  )

  expect_error(
   predict(fit, new_data = new_data[c(), ], pred_horizon = 1000),
   regexp = 'new data are empty'
  )

 }
)



bad_data <- new_data
bad_data$trt <- as.numeric(new_data$trt)

test_that(
 desc = 'unexpected data types are detected',
 code = {
  expect_error(
   object = predict(fit, bad_data, pred_horizon = 1000),
   regexp = "\\<trt\\>"
  )
 }
)

bad_data <- new_data
bad_data$sex <- factor(bad_data$sex, levels = c("m", "f", "new_level"))

test_that(
 desc = 'unexpected factor levels are detected',
 code = {
  expect_error(
   object = predict(fit, bad_data, pred_horizon = 1000),
   regexp = "new_level"
  )
 }
)


bad_data <- new_data
bad_data$sex <- NULL
bad_data$trt <- NULL

test_that(
 desc = 'missing columns are detected',
 code = {
  expect_error(
   object = predict(fit, bad_data, pred_horizon = 1000),
   regexp = "trt and sex"
  )
 }
)

bad_data <- new_data


test_that(
 desc = 'missing values are detected',
 code = {

  bad_data$age[1] <- NA_real_

  expect_error(
   object = predict(fit, bad_data, pred_horizon = 1000),
   regexp = "missing values"
  )

  bad_data$age[1] <- Inf

  expect_error(
   object = predict(fit, bad_data, pred_horizon = 1000),
   regexp = "infinite"
  )


 }
)

test_that(
 desc = 'pred horizon < max time',
 code = {
  expect_error(
   object = predict(fit, pbc_test, pred_horizon = 100000),
   regexp = "max follow-up"
  )
 }
)

test_that(
 desc = "outside limit predictions = predictions at the boundary",
 code = {
  expect_equal(
   predict(fit, pbc_test,
           pred_horizon = 100000,
           boundary_checks = F),
   predict(fit, pbc_test,
           pred_horizon = max(pbc_train$time))
  )
 }
)

test_that(
 desc = 'pred horizon in increasing order',
 code = {

  normal <- predict(fit, pbc_test,
                    pred_horizon = c(2000, 3000, 4000))

  reversed <- predict(fit, pbc_test,
                      pred_horizon = c(4000, 3000, 2000))

  bizaro_1 <- predict(fit, pbc_test,
                      pred_horizon = c(3000, 2000, 4000))

  bizaro_2 <- predict(fit, pbc_test,
                      pred_horizon = c(4000, 2000, 3000))

  bizaro_3 <- predict(fit, pbc_test,
                      pred_horizon = c(3000, 4000, 2000))

  expect_equal(normal, reversed[, c(3,2,1)])
  expect_equal(normal, bizaro_1[, c(2,1,3)])
  expect_equal(normal, bizaro_2[, c(2,3,1)])
  expect_equal(normal, bizaro_3[, c(3,1,2)])

 }
)


test_that(
 desc = 'inconsistent units are detected',
 code = {

  suppressMessages(library(units))
  pbc_units <- pbc_orsf
  units(pbc_units$age) <- 'years'

  pbc_test_units <- pbc_test
  units(pbc_test_units$age) <- 'days'

  fit <- orsf(formula = time + status  ~ . - id,
              data = pbc_units,
              n_tree = n_tree_test)

  expect_error(predict(fit, new_data = pbc_test_units, pred_horizon = 1000),
               'has unit')

 }
)

test_that(
 desc = 'predictions dont require cols in same order as training data',
 code = {

  p1 <- predict(fit, new_data = new_data, pred_horizon = 1000)

  new_col_order <- sample(names(new_data),
                          size = ncol(new_data),
                          replace = F)

  new_data_reordered <- new_data[, new_col_order]

  p2 <- predict(fit, new_data_reordered, pred_horizon = 1000)

  expect_equal(p1, p2)
 }
)


# test_that(
#  'units are vetted in testing data',
#  code = {
#
#   suppressMessages(library(units))
#   pbc_units_trn <- pbc_train
#   pbc_units_tst <- pbc_test
#
#
#   units(pbc_units_trn$time) <- 'days'
#   units(pbc_units_trn$age) <- 'years'
#   units(pbc_units_trn$bili) <- 'mg/dl'
#
#   fit_units = orsf(formula = time + status  ~ . - id,
#                    data = pbc_units_trn,
#                    n_tree = n_tree_test,
#                    oobag_pred_horizon = c(1000, 2500),
#                    tree_seeds = seeds_standard)
#
#   expect_error(
#    predict(fit_units, new_data = pbc_units_tst, pred_horizon = 1000),
#    regexp = 'time, age, and bili'
#   )
#
#   units(pbc_units_tst$time) <- 'years'
#   units(pbc_units_tst$age) <- 'years'
#   units(pbc_units_tst$bili) <- 'mg/dl'
#
#   expect_error(
#    predict(fit_units, new_data = pbc_units_tst, pred_horizon = 1000),
#    regexp = 'time has unit d in the training data'
#   )
#
#   units(pbc_units_tst$time) <- 'days'
#   units(pbc_units_tst$age) <- 'years'
#   units(pbc_units_tst$bili) <- 'mg/dl'
#
#   expect_equal_leaf_summary(fit_units, pred_objects_surv$surv$fit)
#   expect_equal_oobag_eval(fit_units, pred_objects_surv$surv$fit)
#
#   units(pbc_units_tst$time) <- 'days'
#   units(pbc_units_tst$age) <- 'years'
#   units(pbc_units_tst$bili) <- 'mg/l'
#
#   expect_error(
#    predict(fit_units, new_data = pbc_units_tst, pred_horizon = 1000),
#    regexp = 'bili has unit mg/dl in the training data'
#   )
#
#  }
#
# )

# Tests for passing missing data ----

na_index_age <- c(1, 4, 8)
na_index_sex <- c(2, 4, 7)

na_expect <- union(na_index_age, na_index_sex)
obs_expect <- setdiff(1:10, na_expect)

new_data_miss <- pbc_test

new_data_miss$age[na_index_age] <- NA
new_data_miss$sex[na_index_sex] <- NA

new_data_dt_miss <- as.data.table(new_data_miss)
new_data_tbl_miss <- tibble::as_tibble(new_data_miss)

p_cc <- predict(fit,
                new_data = new_data)

p_ps <- predict(fit,
                new_data = new_data_miss,
                na_action = 'pass')

p_ps_dt <- predict(fit,
                   new_data = new_data_dt_miss,
                   na_action = 'pass')

p_ps_tbl <- predict(fit,
                    new_data = new_data_tbl_miss,
                    na_action = 'pass')

test_that(
 desc = "proper error for bad value of na_action",
 code = {
  expect_error(predict(fit,
                       new_data = new_data_miss,
                       na_action = 'failzor'),
               regexp = 'failzor')
 }
)

test_that(
 desc = "same values propagated to pred output with na_action = pass",
 code = {

  expect_equal(p_cc[obs_expect, ],
               p_ps[obs_expect, ],
               tolerance = 0.05)

  expect_equal(p_cc[obs_expect, ],
               p_ps_dt[obs_expect, ],
               tolerance = 0.05)

  expect_equal(p_cc[obs_expect, ],
               p_ps_tbl[obs_expect, ],
               tolerance = 0.05)
 }

)

test_that(
 desc = "missing values propagated to pred output with na_action = pass",
 code = {

  expect_true(all(is.na(p_ps[na_expect, ])))

  expect_equal(p_ps,
               p_ps_dt,
               tolerance = 0.05)

  expect_equal(p_ps,
               p_ps_tbl,
               tolerance = 0.05)
 }
)

# repeat test above with multiple predict horizons

pred_horiz <- c(100, 200, 300, 400, 500)

p_cc <- predict(fit,
                new_data = new_data,
                pred_horizon = pred_horiz)

p_ps <- predict(fit,
                new_data = new_data_miss,
                na_action = 'pass',
                pred_horizon = pred_horiz)

p_ps_dt <- predict(fit,
                   new_data = new_data_dt_miss,
                   na_action = 'pass',
                   pred_horizon = pred_horiz)

p_ps_tbl <- predict(fit,
                    new_data = new_data_tbl_miss,
                    na_action = 'pass',
                    pred_horizon = pred_horiz)

test_that(
 desc = "same values propagated to pred output with na_action = pass",
 code = {

  expect_equal(p_cc[obs_expect, ],
               p_ps[obs_expect, ],
               tolerance = 0.05)

  expect_equal(p_cc[obs_expect, ],
               p_ps_dt[obs_expect, ],
               tolerance = 0.05)

  expect_equal(p_cc[obs_expect, ],
               p_ps_tbl[obs_expect, ],
               tolerance = 0.05)
 }
)

test_that(
 desc = "missing values propagated to pred output with na_action = pass",
 code = {

  expect_true(all(is.na(p_ps[na_expect, ])))

  expect_equal(p_ps,
               p_ps_dt,
               tolerance = 0.05)

  expect_equal(p_ps,
               p_ps_tbl,
               tolerance = 0.05)

 }

)

new_data_all_miss <- new_data_miss

new_data_all_miss$age <- NA_real_

test_that(
 desc = "no blank columns allowed",
 code = {
  expect_error(
   predict(fit, new_data = new_data_all_miss, na_action = 'pass'),
   regexp = 'age has no observed values'
  )
 }
)


