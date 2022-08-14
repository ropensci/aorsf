

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

infer_pred_horizon <- function(object, pred_horizon){

 check_arg_is(object, 'object', 'orsf_fit')

 if(is.null(pred_horizon)) pred_horizon <- object$pred_horizon

 if(is.null(pred_horizon))
  stop("pred_horizon was not specified and could not be found in object. ",
       "did you use oobag_pred = FALSE when running orsf()?",
       call. = FALSE)

 pred_horizon

}
