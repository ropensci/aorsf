

f <- Surv(time, status) ~ . - id


#' @srrstats {G5.2} *Appropriate error behaviour is explicitly demonstrated through tests.*
test_that("inputs are vetted", {

 #' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*

 expect_error(orsf_control_cph(method = 'oh no'), "breslow or efron")

 expect_error(orsf_control_net(alpha = 32), 'should be <= 1')

 f_bad_1 <- function(a_node, y_node, w_node){ 1 }
 f_bad_2 <- function(x_node, a_node, w_node){ 1 }
 f_bad_3 <- function(x_node, y_node, a_node){ 1 }
 f_bad_4 <- function(x_node, y_node){ 1 }

 f_bad_5 <- function(x_node, y_node, w_node) {
  stop("an expected error occurred")
 }

 f_bad_6 <- function(x_node, y_node, w_node){
  return(matrix(0, ncol = 2, nrow = ncol(x_node)))
 }

 f_bad_7 <- function(x_node, y_node, w_node){
  return(matrix(0, ncol = 1, nrow = 2))
 }

 f_bad_8 <- function(x_node, y_node, w_node) {runif(n = ncol(x_node))}

 expect_error(orsf_control_custom(f_bad_1), 'x_node')
 expect_error(orsf_control_custom(f_bad_2), 'y_node')
 expect_error(orsf_control_custom(f_bad_3), 'w_node')
 expect_error(orsf_control_custom(f_bad_4), 'should have 3')
 expect_error(orsf_control_custom(f_bad_5), 'encountered an error')
 expect_error(orsf_control_custom(f_bad_6), 'with 1 column')
 expect_error(orsf_control_custom(f_bad_7), 'with 1 row for each')
 expect_error(orsf_control_custom(f_bad_8), 'matrix output')

 f_rando <- function(x_node, y_node, w_node) { matrix(runif(ncol(x_node)), ncol=1) }

 expect_s3_class(orsf_control_custom(f_rando), 'orsf_control')


})


test_that(
 desc = 'custom orsf_control predictions are good',
 code = {

  fit_pca = orsf(pbc_orsf,
                 Surv(time, status) ~ .,
                 tree_seeds = seeds_standard,
                 control = orsf_control_custom(beta_fun = f_pca),
                 n_tree = n_tree_test)

  expect_gt(fit_pca$eval_oobag$stat_values, .65)

 }
)


