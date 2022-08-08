

object <- orsf(pbc_orsf,
               Surv(time, status) ~ . - id,
               n_tree = 100,
               oobag_pred = TRUE)

n_variables <- 3

smry_1 <- orsf_summarize_uni(object,
                             n_variables = n_variables,
                             pred_type = 'survival')

smry_2 <- orsf_summarize_uni(object,
                             pred_horizon = 1000,
                             n_variables = NULL,
                             pred_type = 'survival')


no_miss_list <- function(l){

 sapply(l, function(x){

  if(is.list(x)) {return(no_miss_list(x))}

  any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x))

 })

}

fi <- get_fctr_info(object)

#' @srrstats {G5.2} *Appropriate error behaviour is explicitly demonstrated through tests.*
#' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*

test_that("output is normal", {


 expect_s3_class(smry_1, class = 'aorsf_summary_uni')
 expect_true(length(unique(smry_1$dt$variable)) == n_variables)
 expect_true(smry_1$pred_horizon == object$pred_horizon)
 expect_true(smry_1$pred_type == 'survival')

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

 expect_error(orsf_summarize_uni(object,
                                 pred_horizon = 1e10,
                                 n_variables = n_variables,
                                 pred_type = 'survival'),
              regexp = 'max time')



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


