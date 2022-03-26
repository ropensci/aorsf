

f <- Surv(time, status) ~ . - id


#' @srrstats {G5.2} *Appropriate error behaviour is explicitly demonstrated through tests.*
test_that("inputs are vetted", {

 #' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*

 expect_error(orsf_control_cph(method = 'oh no'), "breslow or efron")
 expect_error(orsf_control_cph(do_scale = FALSE,
                               iter_max = 10), "must be TRUE")

 expect_error(orsf_control_net(alpha = 32), 'should be <= 1')

})


