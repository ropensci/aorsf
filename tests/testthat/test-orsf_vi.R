

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

oobag_fun_brier <- function(y_mat, s_vec){

 # risk = 1 - survival
 r_vec <- 1 - s_vec

 y <- y_mat[, 'status']

 # mean of the squared differences between predicted and observed risk
 bri <- mean( (y - r_vec)^2 )

 y_mean <- mean(y)

 ref <- mean( (y - y_mean)^2 )

 answer <- 1 - bri / ref

 answer

}

oobag_fun_bad_name <- function(nope, s_vec){

 # risk = 1 - survival
 r_vec <- 1 - s_vec

 # mean of the squared differences between predicted and observed risk
 mean( (y_mat[, 'status'] - r_vec)^2 )

}

oobag_fun_bad_name_2 <- function(y_mat, nope){

 # risk = 1 - survival
 r_vec <- 1 - s_vec

 # mean of the squared differences between predicted and observed risk
 mean( (y_mat[, 'status'] - r_vec)^2 )

}

oobag_fun_bad_out <- function(y_mat, s_vec){

 # risk = 1 - survival
 r_vec <- 1 - s_vec

 # mean of the squared differences between predicted and observed risk
 quantile( (y_mat[, 'status'] - r_vec)^2, probs = c(0.25, 0.50, 0.75) )

}

oobag_fun_bad_out_2 <- function(y_mat, s_vec){

 # mean of the squared differences between predicted and observed risk
 return("A")

}

oobag_fun_3_args <- function(y_mat, s_vec, nope){

 # risk = 1 - survival
 r_vec <- 1 - s_vec

 y <- y_mat[, 'status']

 # mean of the squared differences between predicted and observed risk
 bri <- mean( (y - r_vec)^2 )

 y_mean <- mean(y)

 ref <- mean( (y - y_mean)^2 )

 answer <- 1 - bri / ref

 answer

}


oobag_fun_errors_on_test <- function(y_mat, s_vec){

 stop("I can't do anything!")

}

pbc_vi <- pbc_orsf

pbc_vi$junk <- rnorm(nrow(pbc_orsf))

pbc_vi$junk_cat <- factor(
 sample(letters[1:5], size = nrow(pbc_orsf), replace = TRUE)
)

set.seed(32987)
fit <- orsf(pbc_vi,
            formula = Surv(time, status) ~ age + sex + bili + junk + junk_cat,
            importance = "negate",
            oobag_eval_every = 100)

fit_anova <- orsf(pbc_vi,
                  formula = Surv(time, status) ~ age + sex + bili + junk + junk_cat,
                  importance = "anova",
                  oobag_eval_every = 100)

set.seed(32987)
fit_no_vi <- orsf(pbc_vi,
                  formula = Surv(time, status) ~ age + sex + bili + junk + junk_cat,
                  importance = "none",
                  oobag_eval_every = 100)


test_that(
 desc = 'orsf_vi_ identical to fit importance',
 code = {

  expect_equal(
   fit$importance, # TODO: something going on with name order here.
   orsf_vi_negate(fit, group_factors = FALSE)[names(fit$importance)]
  )

  expect_error(orsf_vi_anova(fit, group_factors = FALSE),
               regexp = "ANOVA")

  expect_equal(
   sort(as.numeric(fit_anova$importance), decreasing = TRUE),
   as.numeric(orsf_vi_anova(fit_anova, group_factors = FALSE))
  )

 }
)

test_that(
 desc = "negation importance picks the right variable",
 code = {
  expect_gt(fit$importance['bili'], fit$importance['junk'])
  expect_gt(fit$importance['bili'], fit$importance['junk_cat_b'])
  expect_gt(fit$importance['bili'], fit$importance['junk_cat_c'])
  expect_gt(fit$importance['bili'], fit$importance['junk_cat_d'])
  expect_gt(fit$importance['bili'], fit$importance['junk_cat_e'])
 }
)

test_that(
 desc = "anova importance picks the right variable",
 code = {
  expect_gt(fit_anova$importance['bili'], fit_anova$importance['junk'])
  expect_gt(fit_anova$importance['bili'], fit_anova$importance['junk_cat_b'])
  expect_gt(fit_anova$importance['bili'], fit_anova$importance['junk_cat_c'])
  expect_gt(fit_anova$importance['bili'], fit_anova$importance['junk_cat_d'])
  expect_gt(fit_anova$importance['bili'], fit_anova$importance['junk_cat_e'])
 }
)

test_that(
 desc = 'cstat from last run of orsf is reproducible',
 code = {

  c_target <- last_value(fit$eval_oobag$stat_values)
  c_estimate <- oobag_c_harrell(
   y_mat = as.matrix(fit$data[, c('time', 'status')]),
   s_vec = fit$surv_oobag
  )

  expect_equal(c_target, c_estimate)

 }
)

test_that(
 desc = 'user defined function returns correct importance values',
 code = {

  expect_equal(
   orsf_vi_negate(fit, group_factors = T),
   orsf_vi_negate(fit, oobag_fun = oobag_c_harrell, group_factors = T)
  )

  expect_equal(
   orsf_vi_negate(fit),
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_c_harrell)
  )

  vi_bri <- orsf_vi_negate(fit, oobag_fun = oobag_fun_brier)

  expect_gt(vi_bri['bili'], vi_bri['junk'])
  expect_gt(vi_bri['bili'], vi_bri['junk_cat'])

 }
)


test_that(
 desc = 'error handling orsf_vi',
 code = {

  expect_error(
   orsf_vi_anova(object = 'nope'),
   regexp = 'inherit'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_bad_name),
   regexp = 'y_mat'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_bad_name_2),
   regexp = 's_vec'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_bad_out),
   regexp = 'length 1'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_bad_out_2),
   regexp = 'type character'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_errors_on_test),
   regexp = 'encountered an error'
  )

  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_3_args),
   regexp = 'has 3'
  )

  fit_no_oob <- orsf(pbc_vi,
                     formula = Surv(time, status) ~ age + sex + bili + junk,
                     oobag_pred = FALSE)

  expect_error(orsf_vi_negate(fit_no_oob), regexp = 'out-of-bag')




 }
)

