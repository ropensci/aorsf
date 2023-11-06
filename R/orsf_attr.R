
# nocov start

#' ORSF attributes
#'
#' these functions are convenient wrappers around attr().
#'
#' @param object an `orsf_fit`
#'
#' @noRd
get_control            <- function(object) attr(object, 'control')
get_mtry               <- function(object) attr(object, 'mtry')


# x attributes
get_rows_x_cc          <- function(object) attr(object, 'rows_x_cc')
get_types_x            <- function(object) attr(object, 'types_x')
get_fctr_info          <- function(object) attr(object, 'fctr_info')
get_na_action          <- function(object) attr(object, 'na_action')
get_numeric_bounds     <- function(object) attr(object, 'numeric_bounds')
get_means              <- function(object) attr(object, 'means')
get_modes              <- function(object) attr(object, 'modes')
get_standard_deviations<- function(object) attr(object, 'standard_deviations')

# y attributes
get_names_y            <- function(object) attr(object, 'names_y')
get_max_time           <- function(object) attr(object, 'max_time')

# shared x & y attributes
get_unit_info          <- function(object) attr(object, 'unit_info')
get_sorted             <- function(object) attr(object, 'sorted')

# counts
get_n_obs              <- function(object) attr(object, 'n_obs')
get_n_tree             <- function(object) attr(object, 'n_tree')
get_n_events           <- function(object) attr(object, 'n_events')
get_n_retry            <- function(object) attr(object, 'n_retry')
get_n_thread           <- function(object) attr(object, 'n_thread')
get_n_split            <- function(object) attr(object, 'n_split')
get_n_leaves_mean      <- function(object) attr(object, 'n_leaves_mean')

# tree node attributes
get_leaf_min_events    <- function(object) attr(object, 'leaf_min_events')
get_leaf_min_obs       <- function(object) attr(object, 'leaf_min_obs')
get_split_min_events   <- function(object) attr(object, 'split_min_events')
get_split_min_obs      <- function(object) attr(object, 'split_min_obs')
get_split_min_stat     <- function(object) attr(object, 'split_min_stat')
get_split_rule         <- function(object) attr(object, 'split_rule')


# oobag attributes
# TODO: global change f_oobag_eval to oobag_eval_fun
get_f_oobag_eval       <- function(object) attr(object, 'f_oobag_eval')
get_oobag_pred         <- function(object) attr(object, 'oobag_pred')
get_oobag_pred_type    <- function(object) attr(object, 'oobag_pred_type')
get_oobag_pred_horizon <- function(object) attr(object, 'oobag_pred_horizon')
get_oobag_eval_every   <- function(object) attr(object, 'oobag_eval_every')
get_oobag_eval_type    <- function(object) attr(object, 'type_oobag_eval')

# importance attributes
get_importance            <- function(object) attr(object, 'importance')
get_importance_values     <- function(object) attr(object, 'importance_values')
get_importance_max_pvalue <- function(object) attr(object, 'importance_max_pvalue')
# TODO: global change group_factors to importance_aggregate
get_group_factors      <- function(object) attr(object, 'group_factors')

# tree attributes
get_tree_type          <- function(object) attr(object, 'tree_type')
get_tree_seeds         <- function(object) attr(object, 'tree_seeds')

# sampler attributes
get_sample_with_replacement <- function(object) attr(object, 'sample_with_replacement')
get_sample_fraction         <- function(object) attr(object, 'sample_fraction')

# misc input attributes
get_weights            <- function(object) attr(object, 'weights_user')
get_verbose_progress   <- function(object) attr(object, 'verbose_progress')


#' ORSF status
#'
#' Determine whether an aorsf model has been trained.
#'
#' @param object an object to check training status for
#'
#' @return logical; `TRUE` if trained, `FALSE` otherwise.
#'
#' @noRd
is_trained <- function(object) return(object$trained)


#' Determine whether object has oobag estimates
#'
#' @param object an object of class 'orsf_fit'
#'
#' @return `TRUE` if oobag survival estimates are present, `FALSE` otherwise
#'
#' @noRd
#'
contains_oobag <- function(object) {!is_empty(object$eval_oobag$stat_values)}

#' Determine whether object has variable importance estimates
#'
#' @param object an object of class 'orsf_fit'
#'
#' @return `TRUE` if variable importance estimates are present, `FALSE` otherwise
#'
#' @noRd
#'
contains_vi <- function(object) {!is_empty(object$importance)}


#' Retrieve x-matrix names
#'
#' @param object an object of class 'orsf_fit'
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

# get the last oobag statistic from an aorsf model
# (last means the most recent)
get_last_oob_stat_value <- function(object){
 object$eval_oobag$stat_values[nrow(object$eval_oobag$stat_values),
                               1,
                               drop = TRUE]
}

# retrieve the lowest ranked variable in terms of importance
get_last_vi <- function(object, name_only = TRUE){
 if(name_only)
  names(last_value(object$importance))
 else
  last_value(object$importance)
}

# nocov end
