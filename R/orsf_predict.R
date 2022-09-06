

#' Compute predictions using ORSF
#'
#' Predicted risk or survival (someday also hazard or mortality)
#'   from an ORSF model.
#'
#' @srrstats {G1.4} *documented with Roxygen*
#' @srrstats {ML1.1} *using the terms 'train' and 'test'.*
#' @srrstats {G2.0a} *specified expectations for length of `pred_horizon`. In general, inputs of length > 1 have the term 'vector' in their description, and inputs of length 1 just have the expected type.*
#' @srrstats {G2.1a} *explicit secondary documentation of expectations on data types of all vector inputs*
 #' @srrstats {G2.8} *As part of initial pre-processing, run checks on inputs to ensure that all other sub-functions receive inputs of a single defined class or type.*
#' @srrstats {ML1.1} *The term 'new_data' are used instead of data_test. There are two reasons for this. First, I am making an effort to be consistent with tidymodels. Second, there is a possibility that users will use predict() without the intention of testing their model, e.g., for interpretation.*
#'
#' @param object (*orsf_fit*) a trained oblique random survival forest
#'   (see [orsf]).
#'
#' @param new_data a `r roxy_data_allowed()` to compute predictions in.
#'
#' @param pred_horizon (_double_) a value or vector indicating the time(s)
#'   that predictions will be calibrated to. E.g., if you were predicting
#'   risk of incident heart failure within the next 10 years, then
#'   `pred_horizon = 10`. `pred_horizon` can be `NULL` if `pred_type` is
#'   `'mort'`, since mortality predictions are aggregated over all
#'   event times
#'
#' @param pred_type (_character_) the type of predictions to compute. Valid
#'   options are
#'
#'   - 'risk' : probability of having an event at or before `pred_horizon`.
#'   - 'surv' : 1 - risk.
#'   - 'chf': cumulative hazard function
#'   - 'mort': mortality prediction
#'
#' @param na_action `r roxy_na_action_header()`
#'
#'   - `r roxy_na_action_fail()`
#'   - `r roxy_na_action_pass()`
#'   - `r roxy_na_action_omit()`
#'
#' @param boundary_checks (_logical_) if `TRUE`, `pred_horizon` will be
#'  checked to make sure the requested values are less than the maximum
#'  observed time in `object`'s training data. If `FALSE`, these checks
#'  are skipped.
#'
#' @param ... `r roxy_dots()`
#'
#' @return a `matrix` of predictions. Column `j` of the matrix corresponds
#'   to value `j` in `pred_horizon`. Row `i` of the matrix corresponds to
#'   row `i` in `new_data`.
#'
#' @details
#'
#' `new_data` must have the same columns with equivalent types as the data
#'   used to train `object`. Also, factors in `new_data` must not have levels
#'   that were not in the data used to train `object`.
#'
#' `pred_horizon` values should not exceed the maximum follow-up time in
#'   `object`'s training data, but if you truly want to do this, set
#'   `boundary_checks = FALSE` and you can use a `pred_horizon` as large
#'   as you want. Note that predictions beyond the maximum follow-up time
#'   in the `object`'s training data are equal to predictions at the
#'   maximum follow-up time, because `aorsf` does not estimate survival
#'   beyond its maximum observed time.
#'
#' If unspecified, `pred_horizon` may be automatically specified as the value
#'   used for `oobag_pred_horizon` when `object` was created (see [orsf]).
#'
#'
#' @export
#'
#' @includeRmd Rmd/orsf_predict_examples.Rmd
#'
predict.orsf_fit <- function(object,
                             new_data,
                             pred_horizon = NULL,
                             pred_type = 'risk',
                             na_action = 'fail',
                             boundary_checks = TRUE,
                             ...){

 # catch any arguments that didn't match and got relegated to ...
 # these arguments are mistaken input names since ... isn't used.
 check_dots(list(...), .f = predict.orsf_fit)

 names_x_data <- intersect(get_names_x(object), names(new_data))

 pred_horizon <- infer_pred_horizon(object, pred_horizon)

 check_predict(object = object,
               new_data = new_data,
               pred_horizon = pred_horizon,
               pred_type = pred_type,
               na_action = na_action,
               boundary_checks = boundary_checks)

 pred_horizon_order <- order(pred_horizon)
 pred_horizon_ordered <- pred_horizon[pred_horizon_order]

 cc <- which(stats::complete.cases(select_cols(new_data, names_x_data)))

 check_complete_cases(cc, na_action, nrow(new_data))

 if(is.null(pred_horizon) && pred_type != 'mort'){
  stop("pred_horizon must be specified for ",
       pred_type, " predictions.", call. = FALSE)
 }

 x_new <- as.matrix(
  ref_code(x_data = new_data[cc, ],
           fi = get_fctr_info(object),
           names_x_data = names_x_data)
 )

 pred_type_cpp <- switch(
  pred_type,
  "risk" = "R",
  "surv" = "S",
  "chf"  = "H",
  "mort" = "M"
 )

 out_values <-
  if(pred_type_cpp == "M"){
   orsf_pred_mort(object, x_new)
  } else if (length(pred_horizon) == 1L) {
   orsf_pred_uni(object$forest, x_new, pred_horizon_ordered, pred_type_cpp)
  } else {
   orsf_pred_multi(object$forest, x_new, pred_horizon_ordered, pred_type_cpp)
  }

 if(na_action == "pass"){

  out <- matrix(nrow = nrow(new_data),
                ncol = length(pred_horizon))

  out[cc, ] <- out_values

 } else {

  out <- out_values

 }

 # output in the same order as pred_horizon
 out[, order(pred_horizon_order), drop = FALSE]

}

orsf_pred_mort <- function(object, x_new){

 pred_mat <- orsf_pred_multi(object$forest,
                             x_new = x_new,
                             time_vec = get_event_times(object),
                             pred_type = 'H')

 matrix(apply(pred_mat, MARGIN = 1, FUN = sum), ncol = 1)

}


