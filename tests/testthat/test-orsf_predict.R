

train <- sample(nrow(pbc_orsf), size = 170)

#' @srrstats {G5.5} *Correctness tests are run with a fixed random seed*
set.seed(1)

aorsf = orsf(formula = time + status  ~ . - id,
             data = pbc_orsf[train, ],
             mtry = 5,
             n_split = 10,
             n_tree = 50,
             oobag_pred = TRUE,
             leaf_min_obs = 15)

new_data <- pbc_orsf[-train, ]

p1 <- predict(aorsf, new_data = new_data, pred_horizon = 1000)
p2 <- predict(aorsf, new_data = new_data, pred_horizon = 2000)

p_multi <- predict(aorsf, new_data = new_data, pred_horizon = c(1000, 2000))

test_that(
 desc = 'multi-time preds are same as uni-time',
 code = {
  expect_equal(cbind(p1, p2), p_multi)
 }
)

test_that(
 desc = 'predictions are bounded',
 code = {
  expect_true(all(p1 <= 1) && all(p1 >= 0))
 })

p2 <- predict(aorsf,
              new_data = new_data,
              pred_horizon = 1000,
              pred_type = 'survival')

test_that(
 desc = 'risk is inverse of survival',
 code = {expect_true(all(p1 == 1 - p2))}
)

#' @srrstats {G5.8, G5.8a} **Edge condition tests** *Zero-length data produce expected behaviour*

test_that(
 desc = 'Boundary case: empty new data throw an error',
 code = {

  expect_error(
   predict(aorsf, new_data = new_data[c(), ], pred_horizon = 1000),
   regexp = 'new data are empty'
  )

  expect_error(
   predict(aorsf, new_data = new_data[c(), ], pred_horizon = 1000),
   regexp = 'new data are empty'
  )

 }
)

#' @srrstats {G5.3} *Test that objects returned contain no missing (`NA`) or undefined (`NaN`, `Inf`) values.*

test_that(
 'No missing, nan, or infinite values in prediction output',
 code = {

  expect_false(any(is.na(p1)))
  expect_false(any(is.nan(p1)))
  expect_false(any(is.infinite(p1)))

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
   object = predict(aorsf, bad_data, pred_horizon = 1000),
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
   object = predict(aorsf, bad_data, pred_horizon = 1000),
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
   object = predict(aorsf, bad_data, pred_horizon = 1000),
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
   object = predict(aorsf, bad_data, pred_horizon = 1000),
   regexp = "missing values"
  )

  bad_data$age[1] <- Inf

  expect_error(
   object = predict(aorsf, bad_data, pred_horizon = 1000),
   regexp = "infinite"
  )


 }
)

test_that(
 desc = 'pred horizon < max time',
 code = {
  expect_error(
   object = predict(aorsf, pbc_orsf[-train,], pred_horizon = 100000),
   regexp = "max follow-up"
  )
 }
)

test_that(
 desc = 'pred horizon in increasing order',
 code = {
  expect_error(
   object = predict(aorsf, pbc_orsf[-train,],
                    pred_horizon = c(4000, 2000)),
   regexp = "ascending"
  )
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

p2 <- predict(aorsf, new_data_reordered, pred_horizon = 1000)

test_that(
 desc = 'predictions dont require cols in same order as training data',
 code = {
  expect_true(
   all(p1 == p2)
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
  set.seed(1)

  fit_units = orsf(formula = time + status  ~ . - id,
                   data = pbc_units_trn,
                   mtry = 5,
                   n_split = 10,
                   n_tree = 50,
                   oobag_pred = TRUE,
                   leaf_min_obs = 15)

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
  units(pbc_units_tst$bili) <- 'mg/l'

  expect_error(
   predict(fit_units, new_data = pbc_units_tst, pred_horizon = 1000),
   regexp = 'bili has unit mg/dl in the training data'
  )

  units(pbc_units_tst$time) <- 'days'
  units(pbc_units_tst$age) <- 'years'
  units(pbc_units_tst$bili) <- 'mg/dl'

  p3 <- predict(fit_units, new_data = pbc_units_tst, pred_horizon = 1000)

  expect_equal(p3, p1)

 }

)




