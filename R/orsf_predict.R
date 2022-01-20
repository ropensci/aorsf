

#' Predict risk or survival
#'
#' @srrstats {G1.4} *documented with Roxygen*
#'
#' @param object (_aorsf_) an oblique random survival forest (ORSF; see [orsf]).
#'
#' @param new_data (_data.frame_) data to compute predictions for. Must have
#'   the same columns with equivalent types as the data used to train `object`.
#'   Also, factors in `new_data` must not have levels that were not in the
#'   data used to train `object`. Last, missing data are not supported.
#'
#' @srrstats {G2.0a} documenting length of `times`.
#'
#' @param times (_double_) a single time indicating the prediction horizon.
#'   Predicted risk or survival values will indicate the probability of
#'   having an event or surviving from baseline to the prediction horizon.
#'   When using [predict.aorsf()], `times` can be a vector of arbitrary length.
#'   When using [orsf_pd_summary()] or [orsf_pd_ice()], `times` must be
#'   length 1. All `times` values must not exceed the maximum follow-up
#'   time in the oblique RSF's training data. Also, `times` must be entered
#'   in ascending order.
#'
#' @param risk (_logical_) if `TRUE`, predicted risk is returned. If `FALSE`,
#'   predicted survival (i.e., 1-risk) is returned.
#'
#'
#' @param ... not used.
#'
#' @return a `matrix` of predictions. Column `j` of the matrix corresponds
#'   to value `j` in `times`. Row `i` of the matrix corresponds to row `i`
#'   in `new_data`.
#'
#' @export
#'
#' @examples
#'
#' train <- seq(1, nrow(pbc_orsf), by = 2)
#' test <- seq(2, nrow(pbc_orsf), by = 2)
#'
#' fit <- orsf(pbc_orsf[train, ], Surv(time, status) ~ . - id)
#'
#' preds <- predict(fit,
#'                  new_data = pbc_orsf[test, ],
#'                  times = c(500, 1500, 2500))
#'
#' head(preds)
#'
#'
predict.aorsf <- function(object,
                          new_data,
                          times,
                          risk = TRUE,
                          ...){

 check_predict(object, new_data, times, risk)

 x_new <- as.matrix(
  ref_code(x_data = new_data,
           fi = get_fctr_info(object),
           names_x_data = get_names_x(object))
 )

 if(length(times) == 1)
  return(orsf_pred_uni(object$forest, x_new, times, risk))

 orsf_pred_multi(object$forest, x_new, times, risk)

}




