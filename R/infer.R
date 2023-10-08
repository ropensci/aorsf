

#' helper for guessing pred_horizon input
#'
#' @param object 'orsf_fit' object
#' @param pred_horizon NULL or a user's specified pred_horizon
#'
#' @return
#'  - if the user gave a pred_horizon, return that.
#'  - else if the object has a pred_horizon, return that
#'  - else throw an error
#'
#' @noRd

infer_pred_horizon <- function(object, pred_type, pred_horizon){

 check_arg_is(object, 'object', 'orsf_fit')

 if(pred_type %in% c("mort", "leaf")){
  # value of pred_horizon does not matter for these types of prediction
  pred_horizon <- 1
 }

 # see if it was previously specified
 if(is.null(pred_horizon)) pred_horizon <- object$pred_horizon

 # throw error if pred_type requires pred_horizon
 if(is.null(pred_horizon)){
  stop("pred_horizon was not specified and could not be found in object.",
       call. = FALSE)
 }


 pred_horizon

}
