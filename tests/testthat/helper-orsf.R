
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

 y <- prep_y(data, names_y_data)
 x <- prep_x(data, fi, names_x_data, means, standard_deviations)
 w <- sample(1:3, nrow(y), replace = TRUE)

 sorted <- collapse::radixorder(y[, 1],  -y[, 2])

 return(
  list(
   x = x[sorted, ],
   y = y[sorted, ],
   w = w[sorted]
  )
 )


}
