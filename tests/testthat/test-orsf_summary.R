

object <- orsf(pbc_orsf,
               Surv(time, status) ~ . - id,
               n_tree = 100,
               oobag_pred = TRUE)

n_variables <- 3
risk <- FALSE

smry <- orsf_summarize_uni(object, n_variables = n_variables, risk = risk)

test_that("output is normal", {
 expect_s3_class(smry, class = 'aorsf_summary_uni')
 expect_true(length(unique(smry$dt$variable)) == n_variables)
 expect_true(smry$times == object$time_pred)
 expect_true(smry$risk == risk)
 expect_invisible(print(smry))
})

test_that("bad inputs caught", {

 expect_error(
  orsf_summarize_uni(object, n_variables = 50, risk = risk),
  "total number of predictors"
 )

})

