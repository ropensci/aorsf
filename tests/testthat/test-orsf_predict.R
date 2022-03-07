

train <- sample(nrow(pbc_orsf), size = 170)

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

p2 <- predict(aorsf, new_data = new_data, pred_horizon = 1000, risk = FALSE)

test_that(
 desc = 'risk is inverse of survival',
 code = {expect_true(all(p1 == 1 - p2))}
)

bad_data <- new_data
bad_data$trt <- as.numeric(new_data$trt)

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
bad_data$age[1] <- NA_real_

test_that(
 desc = 'missing values are detected',
 code = {
  expect_error(
   object = predict(aorsf, bad_data, pred_horizon = 1000),
   regexp = "missing values"
  )
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






