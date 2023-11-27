
#' null operator (copied from rlang)
#' @noRd
`%||%` <-  function (x, y) {
 if (is.null(x))
  y
 else x
}

#' helper for inferring pred_horizon input
#'
#' @param object 'ObliqueForest' object
#' @param pred_horizon NULL or a user's specified pred_horizon
#'
#' @return
#'  - if the user gave a pred_horizon, return that.
#'  - else if the object has a pred_horizon, return that
#'  - else throw an error
#'
#' @noRd

infer_pred_horizon <- function(object, pred_type, pred_horizon){

 check_arg_is(object, 'object', 'ObliqueForest')

 if(pred_type %in% c("mort", "leaf") || object$tree_type != 'survival'){
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


#' helper for inferring outcome type
#'
#' @param names_y_data character vector of outcome names
#' @param data dataset containing outcomes
#'
#' @return character value: 'survival', 'regression' or 'classification'
#'
#' @examples
#'
#' infer_tree_type('bili', pbc_orsf)
#' infer_tree_type('sex', pbc_orsf)
#' infer_tree_type(c('time', 'status'), pbc_orsf)
#' infer_tree_type(Surv(pbc_orsf$time, pbc_orsf$status), pbc_orsf)
#'
#' @noRd
infer_tree_type <- function(names_y_data, data){

 if(length(names_y_data) > 2){
  stop("formula should have at most two variables as the response",
       call. = FALSE)
 }

 if(length(names_y_data) == 2) {
  return("survival")
 }

 if(is.factor(data[[names_y_data]])){
  return("classification")
 } else if(inherits(data[[names_y_data]], 'Surv')) {
  return("survival")
 } else {
  return("regression")
 }

}
