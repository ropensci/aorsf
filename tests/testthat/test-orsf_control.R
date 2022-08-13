

f <- Surv(time, status) ~ . - id


#' @srrstats {G5.2} *Appropriate error behaviour is explicitly demonstrated through tests.*
test_that("inputs are vetted", {

 #' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*

 expect_error(orsf_control_cph(method = 'oh no'), "breslow or efron")

 expect_error(orsf_control_net(alpha = 32), 'should be <= 1')

 f_bad_1 <- function(a_node, y_node, w_node){ 1 }
 f_bad_2 <- function(x_node, a_node, w_node){ 1 }
 f_bad_3 <- function(x_node, y_node, a_node){ 1 }

 expect_error(orsf_control_custom(f_bad_1), 'x_node')
 expect_error(orsf_control_custom(f_bad_2), 'y_node')
 expect_error(orsf_control_custom(f_bad_3), 'w_node')

 f_bad_4 <- function(x_node, y_node, w_node) {runif(n = ncol(x_node))}

 expect_error(orsf_control_custom(f_bad_4), 'matrix output')

 # seems like this one can throw off github actions?
 if (Sys.getenv("run_all_aorsf_tests") == 'yes') {

  f_bad_5 <- function(x_node, y_node, w_node){
   stop("IDK WHAT TO DO", call. = FALSE)
  }

  expect_error(orsf_control_custom(f_bad_5), "encountered an error")

 }




 f <- function(x_node, y_node, w_node) { matrix(runif(ncol(x_node)), ncol=1) }

 expect_s3_class(orsf_control_custom(f), 'aorsf_control')


})


test_that(
 desc = 'outputs meet expectations on prediction accuracy',
 code = {


  f <- function(x_node, y_node, w_node) { matrix(runif(ncol(x_node)), ncol=1) }

  fit_cph = orsf(pbc_orsf,
                  Surv(time, status) ~ .,
                  tree_seeds = seq(500),
                  control = orsf_control_cph(),
                  n_tree = 500)


  fit_rando = orsf(pbc_orsf,
                    Surv(time, status) ~ .,
                    tree_seeds = seq(500),
                    control = orsf_control_custom(beta_fun = f),
                    n_tree = 500)

  expect_lt(fit_rando$eval_oobag$stat_values,
            fit_cph$eval_oobag$stat_values)

  expect_gt(fit_rando$eval_oobag$stat_values, .6)

 }
)


