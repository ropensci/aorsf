

object <- orsf(pbc_orsf,
               Surv(time, status) ~ . - id,
               n_tree = 100)

n_variables <- 3

smry_1 <- orsf_summarize_uni(object,
                             n_variables = n_variables,
                             pred_type = 'surv')


smry_2 <- orsf_summarize_uni(object,
                             pred_horizon = 1000,
                             n_variables = NULL,
                             pred_type = 'surv')

smry_3 <- orsf_summarize_uni(object,
                             pred_type = 'chf')

dt_smry_1 <- as.data.table(smry_1)
dt_smry_2 <- as.data.table(smry_2)
dt_smry_3 <- as.data.table(smry_3)

test_that(
 desc = 'standard summaries run and can be cast to data tables',
 code = {

  expect_true(inherits(dt_smry_1, 'data.table'))
  expect_true(inherits(dt_smry_2, 'data.table'))
  expect_true(inherits(dt_smry_3, 'data.table'))

  expect_gt(nrow(dt_smry_2), nrow(dt_smry_1))

 }


)

no_miss_list <- function(l){

 sapply(l, function(x){

  if(is.list(x)) {return(no_miss_list(x))}

  any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x))

 })

}

fi <- object$get_fctr_info()

test_that("output is normal", {


 expect_s3_class(smry_1, class = 'orsf_summary_uni')
 expect_true(length(unique(smry_1$dt$variable)) == n_variables)
 expect_true(smry_1$pred_horizon == object$pred_horizon)
 expect_true(smry_1$pred_type == 'surv')

 rows_categorical_variables <- smry_1$dt$variable %in% fi$cols
 rows_numeric_variables <- !rows_categorical_variables

 # level should be NA when the variable is numeric
 expect_true(all(is.na(smry_1$dt$level[rows_numeric_variables])))
 # level should not be NA when the variable is categorical
 expect_false(any(is.na(smry_1$dt$level[rows_categorical_variables])))

 rows_categorical_variables <- smry_2$dt$variable %in% fi$cols
 rows_numeric_variables <- !rows_categorical_variables
 # level should be NA when the variable is numeric
 expect_true(all(is.na(smry_2$dt$level[rows_numeric_variables])))
 # level should not be NA when the variable is categorical
 expect_false(any(is.na(smry_2$dt$level[rows_categorical_variables])))

 #' @srrstats {G5.3} *Test that objects returned contain no missing (`NA`) or undefined (`NaN`, `Inf`) values.*
 # only one thing should have missing values (level)
 expect_equal(Reduce(f = sum, x = no_miss_list(smry_1)), 1)
 expect_equal(Reduce(f = sum, x = no_miss_list(smry_2)), 1)


 expect_true(smry_2$pred_horizon == 1000)

})


test_that(
 desc = "print doesn't cause an error",
 code = {
  # we don't need this printing out on the testthat report.
  expect_invisible(p <- capture.output(print(smry_1)))
 }
)



test_that("bad inputs caught", {

 expect_error(
  orsf_summarize_uni(object, n_variables = 50, pred_type = 'risk'),
  "total number of predictors"
 )

})

# high pred horizon
test_that(
 desc = 'higher pred horizon is not allowed for summary',
 code = {

  expect_error(orsf_summarize_uni(fit_standard_pbc$fast,
                                  pred_horizon = 7000),
               regexp = 'prediction horizon')

 }
)

