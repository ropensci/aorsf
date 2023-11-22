
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

# R version written using matrixStats

weighted_variance <- function (x, w = NULL, idxs = NULL,
                               na.rm = FALSE, center = NULL,
                               ...) {
 n <- length(x)
 if (is.null(w)) {
  w <- rep(1, times = n)
 }
 else if (length(w) != n) {
  stop(sprintf("The number of elements in arguments '%s' and '%s' does not match: %.0f != %.0f",
               "w", "x", length(w), n))
 }
 else if (!is.null(idxs)) {
  w <- w[idxs]
 }
 if (!is.null(idxs)) {
  x <- x[idxs]
  n <- length(x)
 }
 na_value <- NA
 storage.mode(na_value) <- storage.mode(x)
 tmp <- (is.na(w) | w > 0)
 if (!all(tmp)) {
  x <- .subset(x, tmp)
  w <- .subset(w, tmp)
  n <- length(x)
 }
 tmp <- NULL
 if (na.rm) {
  keep <- which(!is.na(x))
  x <- .subset(x, keep)
  w <- .subset(w, keep)
  n <- length(x)
  keep <- NULL
 }

 tmp <- is.infinite(w)
 if (any(tmp)) {
  keep <- tmp
  x <- .subset(x, keep)
  n <- length(x)
  w <- rep(1, times = n)
  keep <- NULL
 }
 tmp <- NULL
 if (n <= 1L)
  return(na_value)
 wsum <- sum(w)
 if (is.null(center)) {
  center <- sum(w * x)/wsum
 }
 x <- x - center
 x <- x^2
 lambda <- 1/(wsum - 1)
 sigma2 <- lambda * sum(w * x)
 x <- w <- NULL

 sigma2
}

#' Find cut-point boundaries (R version)
#'
#'  Used to test the cpp version for finding cutpoints
#'
#' @param y_node outcome matrix
#' @param w_node weight vector
#' @param XB linear combination of predictors
#' @param xb_uni unique values in XB
#' @param leaf_min_events min no. of events in a leaf
#' @param leaf_min_obs min no. of observations in a leaf
#'
#' @noRd
#'
#' @return data.frame with description of valid cutpoints
cp_find_bounds_R <- function(y_node,
                             w_node,
                             XB,
                             xb_uni,
                             leaf_min_events,
                             leaf_min_obs){

 status = y_node[, 'status']

 cp_stats <-
  sapply(
   X = xb_uni,
   FUN = function(x){
    c(
     cp = x,
     e_right = sum(status[XB > x] * w_node[XB > x]),
     e_left = sum(status[XB <= x] * w_node[XB <= x]),
     n_right = sum(as.numeric(XB > x) * w_node),
     n_left = sum(as.numeric(XB <= x) * w_node)
    )
   }
  )

 cp_stats <- as.data.frame(t(cp_stats))

 cp_stats$valid_cp = with(
  cp_stats,
   e_right >= leaf_min_events &
   e_left  >= leaf_min_events  &
   n_right >= leaf_min_obs &
   n_left  >= leaf_min_obs
 )

 cp_stats

}

# oobag functions ----

oobag_brier_clsf <- function(y_mat, w_vec, s_vec){
 y_mean <- mean(y_mat)
 1 - mean((y_mat - s_vec)^2) / mean((y_mat - y_mean)^2)
}

oobag_brier_surv <- function(y_mat, w_vec, s_vec){

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

oobag_fun_bad_name_3 <- function(y_mat, nope, s_vec){

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

 col <- sample(ncol(pca$rotation), 1)

 pca$rotation[, col, drop = FALSE]

}

# testing functions ----

expect_equal_leaf_summary <- function(x, y){
 expect_equal(x$forest$leaf_summary,
              y$forest$leaf_summary,
              tolerance = 1e-9)
}

expect_equal_oobag_eval <- function(x, y, tolerance = 1e-9){
 expect_equal(x$eval_oobag$stat_values,
              y$eval_oobag$stat_values,
              tolerance = tolerance)
}

expect_no_missing <- function(x){

 expect_true(!any(is.na(x)))

}

# data processing ----

prep_test_matrices <- function(data, outcomes = c("time", "status")){

 names_y_data <- outcomes
 names_x_data <- setdiff(names(data), outcomes)

 fi <- fctr_info(data, names_x_data)

 types_x_data <- check_var_types(data,
                                 names_x_data,
                                 valid_types = c('numeric',
                                                 'integer',
                                                 'units',
                                                 'factor',
                                                 'ordered'))

 names_x_numeric <- grep(pattern = "^integer$|^numeric$|^units$",
                         x = types_x_data)

 means <- standard_deviations<- modes <- numeric_bounds <- NULL

 numeric_cols <- names_x_data[names_x_numeric]
 nominal_cols <- fi$cols

 if(!is_empty(nominal_cols)){

  modes <- vapply(
   select_cols(data, nominal_cols),
   collapse::fmode,
   FUN.VALUE = integer(1)
  )

 }

 if(!is_empty(numeric_cols)){

  numeric_data <- select_cols(data, numeric_cols)

  numeric_bounds <- matrix(
   data = c(
    collapse::fnth(numeric_data, 0.1),
    collapse::fnth(numeric_data, 0.25),
    collapse::fnth(numeric_data, 0.5),
    collapse::fnth(numeric_data, 0.75),
    collapse::fnth(numeric_data, 0.9)
   ),
   nrow =5,
   byrow = TRUE,
   dimnames = list(c('10%', '25%', '50%', '75%', '90%'),
                   names(numeric_data))
  )

  means <- collapse::fmean(numeric_data)

  standard_deviations <- collapse::fsd(numeric_data)

 }

 if(any(is.na(select_cols(data, names_y_data))))
  stop("Please remove missing values from the outcome variable(s)",
       call. = FALSE)

 cc <- stats::complete.cases(data[, names_x_data])
 data <- data[cc, ]

 if(length(outcomes) > 1){
  y <- prep_y_surv(data, names_y_data)
  sorted <- collapse::radixorder(y[, 1],  -y[, 2])
 } else if(is.factor(data[[names_y_data]])) {
  y <- prep_y_clsf(data, names_y_data)
  sorted <- collapse::seq_row(data)
 } else {
  y <- matrix(data[[names_y_data]], ncol = 1)
  sorted <- collapse::seq_row(data)
 }

 x <- prep_x(data, fi, names_x_data, means, standard_deviations)
 w <- sample(1:3, nrow(y), replace = TRUE)


 return(
  list(
   x = x[sorted, , drop=FALSE],
   y = y[sorted, , drop=FALSE],
   w = w[sorted]
  )
 )


}
