

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

fit_with_vi <- orsf(data_train = pbc_orsf,
                    formula = Surv(time, status) ~ . - id,
                    importance = TRUE,
                    n_tree = 10)

fit_no_vi <- orsf(data_train = pbc_orsf,
                  formula = Surv(time, status) ~ . - id,
                  importance = TRUE,
                  n_tree = 10)

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
        n_tree = 50)
  )

  time_large <- system.time(
   orsf(pbc_orsf,
        control = orsf_control_cph(iter_max = 50, eps = 1e-10),
        Surv(time, status) ~ . -id,
        n_tree = 50)
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
        n_tree = 100)
  )

  expect_true(time_small['elapsed'] < time_large['elapsed'])

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
pbc_temp[, 'bili'] <- NA

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

pbc_temp <- pbc_orsf

noise_vars <- c('bili', 'chol', 'albumin', 'copper', 'alk.phos', 'ast')

pbc_temp[, noise_vars] <- sapply(pbc_temp[, noise_vars], add_noise)

set.seed(89)
fit_orsf <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 10)
set.seed(89)
fit_orsf_2 <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 10)
set.seed(89)
fit_orsf_noised <- orsf(pbc_temp, Surv(time, status) ~ . - id, n_tree = 10)

test_that(
 desc = 'results are identical if a forest is fitted under the same random seed',
 code = {

  # testing a subset of trees for identical betas

  for(i in seq(get_n_tree(fit_orsf))){
   expect_equal(
    object = fit_orsf$forest[[i]]$betas,
    expected = fit_orsf_2$forest[[i]]$betas
   )
  }

 }
)

test_that(
 desc = "results are similar after adding trivial noise",
 code = {

  expect_true(
   max(abs(fit_orsf$signif_means-fit_orsf_noised$signif_means)) < 0.025
  )

  expect_true(
   abs(fit_orsf$eval_oobag$stat_values - fit_orsf_noised$eval_oobag$stat_values) < 0.01
  )

  expect_true(
   max(abs(fit_orsf$surv_oobag-fit_orsf_noised$surv_oobag)) < 0.1
  )

 }

)





