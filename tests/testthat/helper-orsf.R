
library(survival)

# misc functions used for tests ----

no_miss_list <- function(l){

 sapply(l, function(x){

  if(is.list(x)) {return(no_miss_list(x))}

  any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x))

 })

}

add_noise <- function(x, eps = .Machine$double.eps){

 noise <- runif(length(x), min = -eps, max = eps)

 x + noise

}

change_scale <- function(x, mult_by = 1/2){
 x * mult_by
}

# oobag functions ----


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

oobag_fun_errors_on_test <- function(y_mat, w_vec, s_vec){

 stop("expected error occurred!", call. = FALSE)

}

# linear combo functions ----

f_pca <- function(x_node, y_node, w_node) {

 # estimate two principal components.
 pca <- stats::prcomp(x_node, rank. = 2)

 # use a random principal component to split the node
 pca$rotation[, 2, drop = FALSE]

}

# testing functions ----

expect_equal_leaf_summary <- function(x, y){
 expect_equal(x$forest$leaf_summary,
              y$forest$leaf_summary,
              tolerance = 1e-9)
}

expect_equal_oobag_eval <- function(x, y){
 expect_equal(x$eval_oobag$stat_values,
              y$eval_oobag$stat_values,
              tolerance = 1e-9)
}
