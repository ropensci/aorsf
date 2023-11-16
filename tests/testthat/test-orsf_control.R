

f <- Surv(time, status) ~ . - id


#' @srrstats {G5.2} *Appropriate error behaviour is explicitly demonstrated through tests.*
test_that(
 desc = "inputs are vetted",
 code = {

 #' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*

 expect_error(orsf_control_survival(ties = 'oh no'), "breslow or efron")

 expect_error(orsf_control_survival(net_mix = 32), 'should be <= 1')

 f_rando <- function(x_node, y_node, w_node) { matrix(runif(ncol(x_node)), ncol=1) }

 expect_s3_class(orsf_control_survival(method = f_rando), 'orsf_control')


}
)
