

object <- orsf(pbc_orsf,
               Surv(time, status) ~ . - id,
               n_tree = 100,
               oobag_pred = TRUE)

n_variables <- 3
risk <- FALSE

smry_1 <- orsf_summarize_uni(object,
                             n_variables = n_variables,
                             risk = risk)

smry_2 <- orsf_summarize_uni(object,
                             pred_horizon = 1000,
                             n_variables = NULL,
                             risk = risk)

#' @srrstats {G5.2} *Appropriate error behaviour is explicitly demonstrated through tests.*
#' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*

test_that("output is normal", {


 expect_s3_class(smry_1, class = 'aorsf_summary_uni')
 expect_true(length(unique(smry_1$dt$variable)) == n_variables)
 expect_true(smry_1$pred_horizon == object$pred_horizon)
 expect_true(smry_1$risk == risk)



 expect_true(smry_2$pred_horizon == 1000)

 expect_error(orsf_summarize_uni(object,
                                 pred_horizon = 1e10,
                                 n_variables = n_variables,
                                 risk = risk),
              regexp = 'max time')



})

test_that(
 desc = "print doesn't cause an error?",
 code = {
  expect_invisible(print(smry_1))
 }
)

test_that("bad inputs caught", {

 expect_error(
  orsf_summarize_uni(object, n_variables = 50, risk = risk),
  "total number of predictors"
 )

})


