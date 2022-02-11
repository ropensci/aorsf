

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

test_that(
 desc = 'formula inputs are vetted',
 code = {

  expect_error(orsf(pbc_temp, f1), 'not found in data_train')
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
  expect_s3_class(orsf(as.data.table(pbc_orsf), f), 'aorsf')
  expect_s3_class(orsf(as.data.table(pbc_orsf), f,
                       control = orsf_control_net(),
                       n_tree = 1), 'aorsf')
 }
)


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



#
# head(pbc_temp)

