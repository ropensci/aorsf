

#' Predict risk or survival
#'
#' @srrstats {G1.4} *documented with Roxygen*
#'
#' @srrstats {G2.0a} *specified expectations for length of `pred_horizon`. In general, inputs of length > 1 have the term 'vector' in their description, and inputs of length 1 just have the expected type.*
#'
#' @param object (_aorsf_) an oblique random survival forest (ORSF; see [orsf]).
#'
#' @srrstats {ML1.1} *The term 'new_data' are used instead of data_test. There are two reasons for this. First, I am making an effort to be consistent with tidymodels. Second, there is a possibility that users will use predict() without the intention of testing their model, e.g., for interpretation.*
#'
#' @param new_data (_data.frame_) data to compute predictions for. Must have
#'   the same columns with equivalent types as the data used to train `object`.
#'   Also, factors in `new_data` must not have levels that were not in the
#'   data used to train `object`. Last, missing data are not supported.
#'
#' @srrstats {G2.1a} explicit secondary documentation of expectations on data types of all vector inputs
#'
#' @param pred_horizon (_double_) a single time or a vector of times
#'   indicating the time(s) that predicted risk or survival probabilities
#'   will be computed at.
#'
#' @param pred_type (_character_) the type of predictions to return. Valid
#'   options are
#'
#'   - 'risk' : probability of having an event at or before `pred_horizon`.
#'   - 'survival' : 1 - risk.
#'
#' @param ... not used.
#'
#' @return a `matrix` of predictions. Column `j` of the matrix corresponds
#'   to value `j` in `pred_horizon`. Row `i` of the matrix corresponds to
#'   row `i` in `new_data`.
#'
#' @details
#'
#' `pred_horizon` values must not exceed the maximum follow-up time in
#'   `object`'s training data. Also, `pred_horizon` values must be entered
#'   in ascending order.
#'
#' @export
#'
#' @examples
#'
#' #' @srrstats {ML1.1} using the terms 'train' and 'test'.
#'
#' # indices of data used for training the model
#' train <- seq(1, nrow(pbc_orsf), by = 2)
#'
#' # indices of data used to test the trained model.
#' test <- seq(2, nrow(pbc_orsf), by = 2)
#'
#' fit <- orsf(pbc_orsf[train, ], Surv(time, status) ~ . - id)
#'
#' preds <- predict(fit,
#'                  new_data = pbc_orsf[test, ],
#'                  pred_horizon = c(500, 1500, 2500))
#'
#' head(preds)
#'
#'
predict.aorsf <- function(object,
                          new_data,
                          pred_horizon,
                          pred_type = 'risk',
                          ...){

 #' @srrstats {G2.13} *check for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*

 names_x_data <- get_names_x(object)

 if(any(is.na(new_data[, intersect(names_x_data, names(new_data))]))){
  stop("Please remove missing values from new_data, or impute them.",
       call. = FALSE)
 }

 #' @srrstats {G2.8} *As part of initial pre-processing, run checks on inputs to ensure that all other sub-functions receive inputs of a single defined class or type.*

 check_predict(object, new_data, pred_horizon, pred_type)

 x_new <- as.matrix(
  ref_code(x_data = new_data,
           fi = get_fctr_info(object),
           names_x_data = names_x_data)
 )

 risk <- pred_type == 'risk'

 if(length(pred_horizon) == 1L)
  return(orsf_pred_uni(object$forest, x_new, pred_horizon, risk))

 orsf_pred_multi(object$forest, x_new, pred_horizon, risk)

}




