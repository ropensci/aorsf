
library(survival)

# misc functions used for tests

no_miss_list <- function(l){

 sapply(l, function(x){

  if(is.list(x)) {return(no_miss_list(x))}

  any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x))

 })

}

add_noise <- function(x, eps = .Machine$double.eps){

 noise <- rnorm(length(x), mean = 0, sd = eps/2)
 noise <- pmin(noise, eps)
 noise <- pmax(noise, -eps)

 x + noise

}

change_scale <- function(x, mult_by = 1/2){
 x * mult_by
}

f_pca <- function(x_node, y_node, w_node) {

 # estimate two principal components.
 pca <- stats::prcomp(x_node, rank. = 2)

 # use a random principal component to split the node
 pca$rotation[, 2, drop = FALSE]

}

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
