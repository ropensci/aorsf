

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


#' helper for guessing outcome type
#'
#' @param formula formula object
#' @param data data to look for terms in
#'
#' @return character value: 'survival', 'regression' or 'classification'
#'
#' @details
#' formulas without a left hand side will not work. It is assumed
#'   that the formula is checked for an outcome before it gets here.
#'
#'
#' @examples
#'
#' formula_regr <- bili ~ sex + age
#' formula_clsf <- sex ~ bili + age
#' formula_surv <- time + status ~ bili + sex + age
#'
#' infer_outcome_type(formula_regr, pbc_orsf)
#' infer_outcome_type(formula_clsf, pbc_orsf)
#' infer_outcome_type(formula_surv, pbc_orsf)
#'
#' @noRd
infer_outcome_type <- function(formula, data){

 outcome <- as.character(formula[[2]])

 if(length(outcome) >= 2) return("survival")

  if(is.factor(data[[outcome]])){
   return("classification")
  } else if(inherits(data[[outcome]], 'Surv')) {
   return("survival")
  } else {
   return("regression")
  }

  stop("could not infer outcome type", call. = FALSE)

}
