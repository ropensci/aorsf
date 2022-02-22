



#' ORSF attributes
#'
#' these functions are convenient wrappers around attr().
#'
#' @param object an object to get attributes from
#'
#' @return varies depending on which get_ function is used.
#'
#' @noRd
get_max_time        <- function(object) attr(object, 'max_time')
get_n_leaves_mean   <- function(object) attr(object, 'n_leaves_mean')
get_mtry            <- function(object) attr(object, 'mtry')
get_max_time        <- function(object) attr(object, 'max_time')
get_n_obs           <- function(object) attr(object, 'n_obs')
get_n_tree          <- function(object) attr(object, 'n_tree')
get_n_events        <- function(object) attr(object, 'n_events')
get_fctr_info       <- function(object) attr(object, 'fctr_info')
get_unit_info       <- function(object) attr(object, 'unit_info')
get_names_y         <- function(object) attr(object, 'names_y')
get_types_x         <- function(object) attr(object, 'types_x')
get_leaf_min_events <- function(object) attr(object, 'leaf_min_events')
get_leaf_min_obs    <- function(object) attr(object, 'leaf_min_obs')
get_numeric_bounds  <- function(object) attr(object, 'numeric_bounds')


#' Determine whether object has oobag estimates
#'
#' @param object an object of class 'aorsf'
#'
#' @return `TRUE` if oobag survival estimates are present, `FALSE` otherwise
#'
#' @noRd
#'
contains_oobag <- function(object) {!is_empty(object$surv_oobag)}

#' Determine whether object has variable importance estimates
#'
#' @param object an object of class 'aorsf'
#'
#' @return `TRUE` if variable importance estimates are present, `FALSE` otherwise
#'
#' @noRd
#'
contains_vi <- function(object) !is_empty(object$importance)


#' Retrieve x-matrix names
#'
#' @param object an object of class 'aorsf'
#' @param ref_code_names should the names be ref encoded?
#'
#' @return a character vector
#'
#' @noRd

get_names_x <- function(object, ref_code_names = FALSE){
 if(ref_code_names)
  attr(object, "names_x_ref")
 else
  attr(object, 'names_x')
}
