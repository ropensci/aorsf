

object <- fit_standard_pbc$fast

n_variables <- 1

smry <- orsf_summarize_uni(object,
                           n_variables = 1,
                           pred_type = 'surv')



test_that(
 desc = 'summaries can be cast to data tables',
 code = {
  dt_smry <- as.data.table(smry)
  expect_true(inherits(dt_smry, 'data.table'))
 }
)

no_miss_list <- function(l){

 sapply(l, function(x){

  if(is.list(x)) {return(no_miss_list(x))}

  any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x))

 })

}

test_that(desc = "output is normal",
          code = {

 expect_s3_class(smry, class = 'orsf_summary_uni')
 expect_true(length(unique(smry$dt$variable)) == n_variables)
 expect_true(smry$pred_horizon == object$pred_horizon)
 expect_true(smry$pred_type == 'surv')
 expect_equal(Reduce(f = sum, x = no_miss_list(smry)), 0)

})


test_that(
 desc = "print doesn't cause an error",
 code = {
  skip_on_cran()
  # we don't need this printing out on the testthat report.
  expect_invisible(p <- capture.output(print(smry)))
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

  expect_error(orsf_summarize_uni(object,
                                  pred_horizon = 7000),
               regexp = 'prediction horizon')

 }
)

test_that(
 desc = "single class can be specified",
 code = {

  object <- fit_standard_penguin_species$fast

  smry_all <- orsf_summarize_uni(object, n_variables = 1)

  expect_true(all(object$class_levels %in% smry_all$dt$class))

  smry_adelie <- orsf_summarize_uni(object, class = "Adelie", n_variables = 1)

  expect_true(all(smry_adelie$dt$class == "Adelie"))

 }
)
