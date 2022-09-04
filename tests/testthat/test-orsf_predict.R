

#' @srrstats {G5.5} *Correctness tests are run with a fixed random seed*
set.seed(730)

train <- sample(nrow(pbc_orsf), size = 170)

set.seed(329)
fit = orsf(formula = time + status  ~ . - id, data = pbc_orsf[train, ])

set.seed(329)
fit_oobag_risk <- orsf(formula = time + status  ~ . - id,
                       data = pbc_orsf[train, ],
                       oobag_pred_type = 'risk')

set.seed(329)
fit_oobag_chf <- orsf(formula = time + status  ~ . - id,
                      data = pbc_orsf[train, ],
                      oobag_pred_type = 'chf')

test_that(
 desc = 'oobag risk and surv have equivalent C-stat',
 code = {
  expect_equal(
   fit$eval_oobag$stat_values,
   fit_oobag_risk$eval_oobag$stat_values
  )
 }
)



new_data <- pbc_orsf[-train, ]
new_data_dt <- as.data.table(new_data)
new_data_tbl <- tibble::as_tibble(new_data)




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
 desc = 'preds identical with na_action = pass/fail if no missing data',
 code = {
  expect_equal(
   predict(fit, new_data = new_data, na_action = 'fail'),
   predict(fit, new_data = new_data, na_action = 'pass')
  )
 }
)

p1_risk <- predict(fit, new_data = new_data, pred_horizon = 1000)

p2_risk <- predict(fit_oobag_risk, new_data = new_data, pred_horizon = 1000)

p3_risk <- predict(fit_oobag_chf, new_data = new_data, pred_horizon = 1000)

test_that(
 desc = "same predictions from the forest regardless of oob type",
 code = {
  expect_equal(p1_risk, p2_risk)
  expect_equal(p1_risk, p3_risk)
 }
)

p1_chf <- predict(fit, new_data = new_data,
                  pred_type = 'chf', pred_horizon = 1000)

p1_surv <- predict(fit, new_data = new_data,
                   pred_type = 'surv', pred_horizon = 1000)

p1_mort <- predict(fit, new_data = new_data, pred_type = 'mort')

test_that(
 desc = 'predict is type stable',
 code = {
  expect_equal(dim(p1_risk), dim(p1_chf))
  expect_equal(dim(p1_risk), dim(p1_surv))
  expect_equal(dim(p1_risk), dim(p1_mort))
  expect_equal(dim(p1_chf), dim(p1_surv))
  expect_equal(dim(p1_chf), dim(p1_mort))
  expect_equal(dim(p1_surv), dim(p1_mort))
 }
)


test_that(
 desc = 'predictions computed for tibbles, and data.tables',
 code = {

  p1_dt <- predict(fit, new_data = new_data_dt, pred_horizon = 1000)
  p1_tbl <- predict(fit, new_data = new_data_tbl, pred_horizon = 1000)

  expect_equal(p1_risk, p1_dt)
  expect_equal(p1_risk, p1_tbl)

  p1_dt <- predict(fit, new_data = new_data_dt, pred_type = 'mort')
  p1_tbl <- predict(fit, new_data = new_data_tbl, pred_type = 'mort')

  expect_equal(p1_mort, p1_dt)
  expect_equal(p1_mort, p1_tbl)

  p1_dt <- predict(fit, new_data = new_data_dt,
                   pred_type = 'chf', pred_horizon = 1000)

  p1_tbl <- predict(fit, new_data = new_data_tbl,
                    pred_type = 'chf', pred_horizon = 1000)

  expect_equal(p1_chf, p1_dt)
  expect_equal(p1_chf, p1_tbl)

  p1_dt <- predict(fit, new_data = new_data_dt,
                   pred_type = 'surv', pred_horizon = 1000)

  p1_tbl <- predict(fit, new_data = new_data_tbl,
                    pred_type = 'surv', pred_horizon = 1000)

  expect_equal(p1_surv, p1_dt)
  expect_equal(p1_surv, p1_tbl)


 }

)


p2 <- predict(fit, new_data = new_data, pred_horizon = 2000)

p_multi <- predict(fit, new_data = new_data, pred_horizon = c(1000, 2000))

test_that(
 desc = 'multi-time preds are same as uni-time',
 code = {
  expect_equal(cbind(p1_risk, p2), p_multi)
 }
)

test_that(
 desc = 'predictions are bounded',
 code = {
  expect_true(all(p1_risk <= 1) && all(p1_risk >= 0))
 })

p2 <- predict(fit,
              new_data = new_data,
              pred_horizon = 1000,
              pred_type = 'surv')

test_that(
 desc = 'risk is inverse of survival',
 code = {expect_true(all(p1_risk == 1 - p2))}
)

test_that(
 desc = 'predictions do not depend on observations in the data',
 code = {

  for(i in seq(nrow(new_data))){
   p2_1row <- predict(fit,
                      new_data = new_data[i,],
                      pred_horizon = 1000,
                      pred_type = 'surv')

   expect_equal(p2_1row, p2[i], ignore_attr = TRUE, tolerance = 0.015)
  }
 }
)


set.seed(329)

test_that(
 'predictions do not depend on order of the data',
 code = {

  new_order <- sample(nrow(new_data), replace = F)

  preds <- predict(fit,
                   new_data = new_data[new_order, ],
                   pred_horizon = 1000,
                   pred_type = 'surv')

  expect_equal(preds, p2[new_order], ignore_attr = TRUE, tolerance = 0.015)

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

#' @srrstats {G5.8, G5.8a} **Edge condition tests** *Zero-length data produce expected behaviour*

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

#' @srrstats {G5.3} *Test that objects returned contain no missing (`NA`) or undefined (`NaN`, `Inf`) values.*

test_that(
 'No missing, nan, or infinite values in prediction output',
 code = {

  expect_false(any(is.na(p1_risk)))
  expect_false(any(is.nan(p1_risk)))
  expect_false(any(is.infinite(p1_risk)))

  expect_false(any(is.na(p2)))
  expect_false(any(is.nan(p2)))
  expect_false(any(is.infinite(p2)))

  expect_false(any(is.na(p_multi)))
  expect_false(any(is.nan(p_multi)))
  expect_false(any(is.infinite(p_multi)))

 }
)

bad_data <- new_data
bad_data$trt <- as.numeric(new_data$trt)

#' @srrstats {G5.2} *Appropriate error behaviour is explicitly demonstrated through tests.*
#' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*

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
   object = predict(fit, pbc_orsf[-train,], pred_horizon = 100000),
   regexp = "max follow-up"
  )
 }
)

test_that(
 desc = "outside limit predictions = predictions at the boundary",
 code = {
  expect_equal(
   predict(fit, pbc_orsf[-train,],
           pred_horizon = 100000,
           boundary_checks = F),
   predict(fit, pbc_orsf[-train,],
           pred_horizon = get_max_time(fit))
  )
 }
)

test_that(
 desc = 'pred horizon in increasing order',
 code = {

  normal <- predict(fit, pbc_orsf[-train,],
                    pred_horizon = c(2000, 3000, 4000))

  reversed <- predict(fit, pbc_orsf[-train,],
                      pred_horizon = c(4000, 3000, 2000))

  bizaro_1 <- predict(fit, pbc_orsf[-train,],
                      pred_horizon = c(3000, 2000, 4000))

  bizaro_2 <- predict(fit, pbc_orsf[-train,],
                      pred_horizon = c(4000, 2000, 3000))

  bizaro_3 <- predict(fit, pbc_orsf[-train,],
                      pred_horizon = c(3000, 4000, 2000))

  expect_equal(normal, reversed[, c(3,2,1)])
  expect_equal(normal, bizaro_1[, c(2,1,3)])
  expect_equal(normal, bizaro_2[, c(2,3,1)])
  expect_equal(normal, bizaro_3[, c(3,1,2)])

 }
)




test_that(
 desc = 'missing units are detected',
 code = {

  suppressMessages(library(units))
  pbc_units <- pbc_orsf
  units(pbc_units$age) <- 'years'

  fit <- orsf(formula = time + status  ~ . - id,
              data = pbc_units,
              n_tree = 10)

  expect_error(predict(fit, new_data = pbc_orsf, pred_horizon = 1000),
               'unit attributes')

 }
)

new_col_order <- sample(names(new_data),
                        size = ncol(new_data),
                        replace = F)

new_data_reordered <- new_data[, new_col_order]

p2 <- predict(fit, new_data_reordered, pred_horizon = 1000)

test_that(
 desc = 'predictions dont require cols in same order as training data',
 code = {
  expect_true(
   all(p1_risk == p2)
  )
 }
)


#' @srrstats {G2.11} *test to make sure testing units are consistent with training units when someone is trying to compute predictions.*

test_that(
 'units are vetted in testing data',
 code = {

  suppressMessages(library(units))
  pbc_units_trn <- pbc_orsf[train, ]
  pbc_units_tst <- pbc_orsf[-train, ]


  units(pbc_units_trn$time) <- 'days'
  units(pbc_units_trn$age) <- 'years'
  units(pbc_units_trn$bili) <- 'mg/dl'

  #' @srrstats {G5.5} *Correctness tests are run with a fixed random seed*

  set.seed(329)
  fit_units = orsf(formula = time + status  ~ . - id, data = pbc_units_trn)

  expect_error(
   predict(fit_units, new_data = pbc_units_tst, pred_horizon = 1000),
   regexp = 'time, age, and bili'
  )

  units(pbc_units_tst$time) <- 'years'
  units(pbc_units_tst$age) <- 'years'
  units(pbc_units_tst$bili) <- 'mg/dl'

  expect_error(
   predict(fit_units, new_data = pbc_units_tst, pred_horizon = 1000),
   regexp = 'time has unit d in the training data'
  )

  units(pbc_units_tst$time) <- 'days'
  units(pbc_units_tst$age) <- 'years'
  units(pbc_units_tst$bili) <- 'mg/dl'


  expect_equal(fit_units$pred_oobag, fit$pred_oobag)
  expect_equal(fit_units$eval_oobag$stat_values, fit$eval_oobag$stat_values)
  expect_equal(fit_units$forest, fit$forest)

  # small difference in one or two cases, but the forests are identical.
  # so...
  p3 <- predict(fit_units, new_data = pbc_units_tst, pred_horizon = 1000)
  expect_equal(p3, p1_risk)

  units(pbc_units_tst$time) <- 'days'
  units(pbc_units_tst$age) <- 'years'
  units(pbc_units_tst$bili) <- 'mg/l'

  expect_error(
   predict(fit_units, new_data = pbc_units_tst, pred_horizon = 1000),
   regexp = 'bili has unit mg/dl in the training data'
  )

 }

)

# Tests for passing missing data ----

na_index_age <- c(1, 4, 8)
na_index_sex <- c(2, 4, 7)

na_expect <- union(na_index_age, na_index_sex)
obs_expect <- setdiff(1:10, na_expect)

new_data_miss <- pbc_orsf[-train, ]

new_data_miss$age[na_index_age] <- NA
new_data_miss$sex[na_index_sex] <- NA

new_data_dt_miss <- as.data.table(new_data_miss)
new_data_tbl_miss <- tibble::as_tibble(new_data_miss)

p_cc <- predict(fit,
                new_data = new_data[1:10, ])

p_ps <- predict(fit,
                new_data = new_data_miss[1:10, ],
                na_action = 'pass')

p_ps_dt <- predict(fit,
                   new_data = new_data_dt_miss[1:10, ],
                   na_action = 'pass')

p_ps_tbl <- predict(fit,
                    new_data = new_data_tbl_miss[1:10, ],
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
  expect_identical(p_cc[obs_expect, ],
                   p_ps[obs_expect, ])

  expect_identical(p_cc[obs_expect, ],
                   p_ps_dt[obs_expect, ])

  expect_identical(p_cc[obs_expect, ],
                   p_ps_tbl[obs_expect, ])
 }
)

test_that(
 desc = "missing values propagated to pred output with na_action = pass",
 code = {
  expect_true(all(is.na(p_ps[na_expect, ])))
  expect_identical(p_ps, p_ps_dt)
  expect_identical(p_ps, p_ps_tbl)
 }
)

# repeat test above with multiple predict horizons

pred_horiz <- c(100, 200, 300, 400, 500)

p_cc <- predict(fit,
                new_data = new_data[1:10, ],
                pred_horizon = pred_horiz)

p_ps <- predict(fit,
                new_data = new_data_miss[1:10, ],
                na_action = 'pass',
                pred_horizon = pred_horiz)

p_ps_dt <- predict(fit,
                   new_data = new_data_dt_miss[1:10, ],
                   na_action = 'pass',
                   pred_horizon = pred_horiz)

p_ps_tbl <- predict(fit,
                    new_data = new_data_tbl_miss[1:10, ],
                    na_action = 'pass',
                    pred_horizon = pred_horiz)

test_that(
 desc = "same values propagated to pred output with na_action = pass",
 code = {
  expect_identical(p_cc[obs_expect, ],
                   p_ps[obs_expect, ])

  expect_identical(p_cc[obs_expect, ],
                   p_ps_dt[obs_expect, ])

  expect_identical(p_cc[obs_expect, ],
                   p_ps_tbl[obs_expect, ])
 }
)

test_that(
 desc = "missing values propagated to pred output with na_action = pass",
 code = {
  expect_true(all(is.na(p_ps[na_expect, ])))
  expect_identical(p_ps, p_ps_dt)
  expect_identical(p_ps, p_ps_tbl)
 }
)

new_data_all_miss <- new_data_miss

new_data_all_miss$age <- NA_real_

test_that(
 desc = "can't give orsf nothing but missing data",
 code = {
  expect_error(
   predict(fit, new_data = new_data_all_miss, na_action = 'pass'),
   regexp = 'complete data'
  )
 }
)


