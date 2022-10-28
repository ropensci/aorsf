



#' ORSF attributes
#'
#' these functions are convenient wrappers around attr().
#'
#' @param object an object to get attributes from
#'
#' @return varies depending on which get_ function is used.
#'
#' @noRd
get_mtry               <- function(object) attr(object, 'mtry')
get_n_obs              <- function(object) attr(object, 'n_obs')
get_n_tree             <- function(object) attr(object, 'n_tree')
get_names_y            <- function(object) attr(object, 'names_y')
get_types_x            <- function(object) attr(object, 'types_x')
get_n_events           <- function(object) attr(object, 'n_events')
get_max_time           <- function(object) attr(object, 'max_time')
get_unit_info          <- function(object) attr(object, 'unit_info')
get_fctr_info          <- function(object) attr(object, 'fctr_info')
get_n_leaves_mean      <- function(object) attr(object, 'n_leaves_mean')
get_n_split            <- function(object) attr(object, 'n_split')
get_leaf_min_events    <- function(object) attr(object, 'leaf_min_events')
get_leaf_min_obs       <- function(object) attr(object, 'leaf_min_obs')
get_split_min_events   <- function(object) attr(object, 'split_min_events')
get_split_min_obs      <- function(object) attr(object, 'split_min_obs')
get_split_min_stat     <- function(object) attr(object, 'split_min_stat')
get_cph_method         <- function(object) attr(object, 'cph_method')
get_cph_eps            <- function(object) attr(object, 'cph_eps')
get_cph_iter_max       <- function(object) attr(object, 'cph_iter_max')
get_cph_do_scale       <- function(object) attr(object, 'cph_do_scale')
get_net_alpha          <- function(object) attr(object, 'net_alpha')
get_net_df_target      <- function(object) attr(object, 'net_df_target')
get_numeric_bounds     <- function(object) attr(object, 'numeric_bounds')
get_n_retry            <- function(object) attr(object, 'n_retry')
get_f_oobag_eval       <- function(object) attr(object, 'f_oobag_eval')
get_type_oobag_eval    <- function(object) attr(object, 'type_oobag_eval')
get_oobag_pred         <- function(object) attr(object, 'oobag_pred')
get_oobag_pred_type    <- function(object) attr(object, 'oobag_pred_type')
get_oobag_pred_horizon <- function(object) attr(object, 'oobag_pred_horizon')
get_oobag_eval_every   <- function(object) attr(object, 'oobag_eval_every')
get_importance         <- function(object) attr(object, 'importance')
get_f_beta             <- function(object) attr(object, 'f_beta')
get_orsf_type          <- function(object) attr(object, 'orsf_type')
get_f_oobag_eval       <- function(object) attr(object, 'f_oobag_eval')
get_type_oobag_eval    <- function(object) attr(object, 'type_oobag_eval')
get_tree_seeds         <- function(object) attr(object, 'tree_seeds')
get_weights_user       <- function(object) attr(object, 'weights_user')
get_event_times        <- function(object) attr(object, 'event_times')
get_impute_values      <- function(object) attr(object, 'impute_values')
get_standard_deviations<- function(object) attr(object, 'standard_deviations')
#' ORSF status
#'
#' Determine whether an aorsf model has been trained.
#'
#' @param object an object to check training status for
#'
#' @return logical; `TRUE` if trained, `FALSE` otherwise.
#'
#' @noRd
is_trained <- function(object) attr(object, 'trained')


#' Determine whether object has oobag estimates
#'
#' @param object an object of class 'orsf_fit'
#'
#' @return `TRUE` if oobag survival estimates are present, `FALSE` otherwise
#'
#' @noRd
#'
contains_oobag <- function(object) {!is_empty(object$pred_oobag)}

#' Determine whether object has variable importance estimates
#'
#' @param object an object of class 'orsf_fit'
#'
#' @return `TRUE` if variable importance estimates are present, `FALSE` otherwise
#'
#' @noRd
#'
contains_vi <- function(object) get_importance(object) != 'none'


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
