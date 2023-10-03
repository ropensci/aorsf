
require(survival)

# set id to a factor so that it can trigger the id error
pbc_orsf$id <- factor(pbc_orsf$id)
pbc_orsf$status <- pbc_orsf$status+1

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

  expect_error(orsf(pbc_orsf, f1), 'not found in data')
  # # dropped - see https://github.com/mlr-org/mlr3extralearners/issues/259
  # expect_warning(orsf(pbc_orsf, f2), 'at least 2 predictors')
  expect_error(orsf(pbc_orsf, f3), 'unrecognized')
  expect_error(orsf(pbc_orsf, f4), 'unrecognized')
  expect_error(orsf(pbc_orsf, f5), 'id variable?')
  expect_error(orsf(pbc_orsf, f6), 'not_right')
  expect_error(orsf(pbc_orsf, f7), 'not_right')
  expect_error(orsf(pbc_orsf, f8), 'must have two variables')
  expect_error(orsf(pbc_orsf, f9), 'Did you enter')
  expect_error(orsf(pbc_orsf, f10), 'must have two variables')
  expect_error(orsf(pbc_orsf, f11), 'detected >1 event type')
  expect_error(orsf(pbc_orsf, f13), 'must be two sided')
  expect_error(orsf(pbc_orsf, f14), 'Did you enter')
  expect_error(orsf(pbc_orsf, f15), "as many levels as there are rows")

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

  fit_long <- orsf(formula = f_long, pbc_orsf, n_tree = 10)

  # fits the orsf as expected
  expect_s3_class(fit_long, 'orsf_fit')
  # keeps unique names
  expect_equal(x_vars, get_names_x(fit_long))

 }
)
