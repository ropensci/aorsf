

oobag_fun_brier <- function(y_mat, w_vec, s_vec){

 # risk = 1 - survival
 r_vec <- 1 - s_vec

 y <- y_mat[, 2L]

 # mean of the squared differences between predicted and observed risk
 bri <- mean( (y - r_vec)^2 )

 y_mean <- mean(y)

 ref <- mean( (y - y_mean)^2 )

 answer <- 1 - bri / ref

 answer

}

oobag_fun_bad_name <- function(nope, w_vec, s_vec){

 # risk = 1 - survival
 r_vec <- 1 - s_vec

 # mean of the squared differences between predicted and observed risk
 mean( (y_mat[, 2L] - r_vec)^2 )

}

oobag_fun_bad_name_2 <- function(y_mat, w_vec, nope){

 # risk = 1 - survival
 r_vec <- 1 - s_vec

 # mean of the squared differences between predicted and observed risk
 mean( (y_mat[, 2L] - r_vec)^2 )

}

oobag_fun_bad_out <- function(y_mat, w_vec, s_vec){

 # risk = 1 - survival
 r_vec <- 1 - s_vec

 # mean of the squared differences between predicted and observed risk
 quantile( (y_mat[, 2L] - r_vec)^2, probs = c(0.25, 0.50, 0.75) )

}

oobag_fun_bad_out_2 <- function(y_mat, w_vec, s_vec){

 # mean of the squared differences between predicted and observed risk
 return("A")

}

oobag_fun_4_args <- function(y_mat, w_vec, s_vec, nope){

 # risk = 1 - survival
 r_vec <- 1 - s_vec

 y <- y_mat[, 2L]

 # mean of the squared differences between predicted and observed risk
 bri <- mean( (y - r_vec)^2 )

 y_mean <- mean(y)

 ref <- mean( (y - y_mean)^2 )

 answer <- 1 - bri / ref

 answer

}


oobag_fun_errors_on_test <- function(y_mat, s_vec){

 stop("I can't do anything!", call. = FALSE)

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
            group_factors = FALSE,
            oobag_eval_every = 100)

set.seed(32987)
fit_anova <- orsf(pbc_vi,
                  formula = Surv(time, status) ~ age + sex + bili + junk + junk_cat,
                  importance = "anova",
                  group_factors = FALSE,
                  oobag_eval_every = 100)

set.seed(32987)
fit_permute <- orsf(pbc_vi,
                    formula = Surv(time, status) ~ age + sex + bili + junk + junk_cat,
                    importance = "permute",
                    group_factors = FALSE,
                    oobag_eval_every = 100)

set.seed(32987)
fit_no_vi <- orsf(pbc_vi,
                  formula = Surv(time, status) ~ age + sex + bili + junk + junk_cat,
                  importance = "none",
                  group_factors = FALSE,
                  oobag_eval_every = 100)


test_that(
 desc = 'orsf_vi_ identical to fit importance',
 code = {

  expect_equal(
   fit$importance, # something going on with name order here.
   orsf_vi_negate(fit, group_factors = FALSE)[names(fit$importance)]
  )

  expect_equal(
   fit$importance,
   orsf_vi(fit,
           group_factors = FALSE,
           importance = 'negate')[names(fit$importance)]
  )

  # can't extract what isn't there
  expect_error(orsf_vi(fit_no_vi), regexp = 'no variable importance')

  # general api function matches the expected values of lower-level things
  expect_equal(orsf_vi(fit_anova, group_factors = F), fit_anova$importance)

  # negation results identical across api funs
  expect_equal(
   orsf_vi_negate(fit_no_vi),
   orsf_vi(fit_no_vi, importance = 'negate')
  )

  # permutation results identical across api funs using same seed
  set.seed(329)
  vi_permute_1 <- orsf_vi_permute(fit_no_vi)
  set.seed(329)
  vi_permute_2 <- orsf_vi(fit_no_vi, importance = 'permute')

  expect_equal(vi_permute_2, vi_permute_1)

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
 desc = "permutation importance picks the right variable",
 code = {
  expect_gt(fit_permute$importance['bili'], fit_anova$importance['junk'])
  expect_gt(fit_permute$importance['bili'], fit_anova$importance['junk_cat_b'])
  expect_gt(fit_permute$importance['bili'], fit_anova$importance['junk_cat_c'])
  expect_gt(fit_permute$importance['bili'], fit_anova$importance['junk_cat_d'])
  expect_gt(fit_permute$importance['bili'], fit_anova$importance['junk_cat_e'])
 }
)

test_that(
 desc = 'cstat from last run of orsf is reproducible',
 code = {

  c_target <- last_value(fit$eval_oobag$stat_values)
  c_estimate <- oobag_c_survival(
   y_mat = as.matrix(fit$data[, c('time', 'status')]),
   w_vec = rep(1, nrow(fit$data)),
   s_vec = fit$pred_oobag
  )

  expect_equal(c_target, c_estimate)

 }
)

test_that(
 desc = 'user defined function returns correct importance values',
 code = {

  expect_equal(
   orsf_vi_negate(fit, group_factors = T),
   orsf_vi_negate(fit, oobag_fun = oobag_c_risk, group_factors = T)
  )

  expect_equal(
   orsf_vi_negate(fit),
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_c_risk)
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

  if(Sys.getenv("run_all_aorsf_tests") == 'yes'){
   expect_error(
    orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_errors_on_test),
    regexp = 'encountered an error'
   )
  }


  expect_error(
   orsf_vi_negate(fit_no_vi, oobag_fun = oobag_fun_4_args),
   regexp = 'has 4'
  )

 }
)

data_with_empty_factor <- pbc_orsf

levels(data_with_empty_factor$sex) <-
 c(levels(data_with_empty_factor$sex), 'o')

fit <- orsf(data_with_empty_factor,
            formula = time + status ~ . - id,
            group_factors = FALSE)

test_that(
 desc = 'unused factor levels have nan for importance, but this gets kicked out of the aggregated factor importance',
 code = {
  expect_equal(sum(is.nan(fit$importance)), 1)
  expect_equal(sum(is.nan(orsf_vi(fit))), 0)
 }
)


# getting ungrouped VI values from an orsf fit


test_that(
 desc = "ungrouped VI can be recovered with group_factors = TRUE in orsf()",
 code = {

  expect_true('sex' %in% names(orsf_vi(fit_anova)))
  expect_true('sex_f' %in% names(orsf_vi(fit_anova, group_factors = FALSE)))

  permute_importance <- orsf_vi_permute(fit_anova, group_factors = FALSE)
  negate_importance <- orsf_vi_negate(fit_anova, group_factors = FALSE)

  expect_true('sex_f' %in% names(permute_importance))
  expect_true('sex_f' %in% names(negate_importance))

  # same but with group_factors TRUE

  permute_importance <- orsf_vi_permute(fit_anova, group_factors = TRUE)
  negate_importance <- orsf_vi_negate(fit_anova, group_factors = TRUE)

  expect_true('sex' %in% names(permute_importance))
  expect_true('sex' %in% names(negate_importance))


 }
)

# repeat some tests with grouped factors
set.seed(32987)
fit <- orsf(pbc_vi,
            formula = Surv(time, status) ~ age + sex + bili + junk + junk_cat,
            importance = "negate",
            group_factors = TRUE,
            oobag_eval_every = 100)

set.seed(32987)
fit_anova <- orsf(pbc_vi,
                  formula = Surv(time, status) ~ age + sex + bili + junk + junk_cat,
                  importance = "anova",
                  group_factors = TRUE,
                  oobag_eval_every = 100)

set.seed(32987)
fit_permute <- orsf(pbc_vi,
                    formula = Surv(time, status) ~ age + sex + bili + junk + junk_cat,
                    importance = "permute",
                    group_factors = TRUE,
                    oobag_eval_every = 100)

set.seed(32987)
fit_no_vi <- orsf(pbc_vi,
                  formula = Surv(time, status) ~ age + sex + bili + junk + junk_cat,
                  importance = "none",
                  group_factors = TRUE,
                  oobag_eval_every = 100)

test_that(
 desc = "negation importance picks the right variable",
 code = {
  expect_gt(fit$importance['bili'], fit$importance['junk'])
  expect_gt(fit$importance['bili'], fit$importance['junk_cat'])
 }
)

test_that(
 desc = "anova importance picks the right variable",
 code = {
  expect_gt(fit_anova$importance['bili'], fit_anova$importance['junk'])
  expect_gt(fit_anova$importance['bili'], fit_anova$importance['junk_cat'])
 }
)

test_that(
 desc = "permutation importance picks the right variable",
 code = {
  expect_gt(fit_permute$importance['bili'], fit_anova$importance['junk'])
  expect_gt(fit_permute$importance['bili'], fit_anova$importance['junk_cat'])
 }
)

test_that(
 desc = "can get both types of importance from scratch, grouped or not",
 code = {

  permute <- orsf_vi_permute(fit_no_vi, group_factors = FALSE)
  negate <- orsf_vi_negate(fit_no_vi, group_factors = FALSE)

  expect_true('sex_f' %in% names(permute))
  expect_true('sex_f' %in% names(negate))

  permute <- orsf_vi_permute(fit_no_vi, group_factors = TRUE)
  negate <- orsf_vi_negate(fit_no_vi, group_factors = TRUE)

  expect_true('sex' %in% names(permute))
  expect_true('sex' %in% names(negate))


 }
)

