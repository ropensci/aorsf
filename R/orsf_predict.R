

#' Prediction with oblique RSF
#'
#' @param object (_aorsf_) an oblique random survival forest (RSF; see [orsf]).
#'
#' @param new_data (_data.frame_) data to compute predictions for. Must have
#'   the same columns with equivalent types as the data used to train `object`.
#'   Also, factors in `new_data` must not have levels that were not in the
#'   data used to train `object`. Last, missing data are not supported.
#'
#' @param times (_double_) a single time or a vector of times for oblique RSF
#'   predictions. All `times` values must not exceed the maximum follow-up
#'   time in the oblique RSF's training data. Also, `times` must be entered
#'   in ascending order.
#'
#' @param risk (_logical_) if `TRUE`, predicted risk is returned. If `FALSE`,
#'   predicted survival (i.e., 1-risk) is returned.
#'
#' @param run_checks (_logical_) If `TRUE`, tests will be run on inputs.
#'  If `FALSE`, the bare minimum will be checked. It is very easy
#'  to crash your R session when `run_checks` is `FALSE`. Use with care!
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
predict.aorsf <- function(object, new_data, times,
                          risk = TRUE,
                          run_checks = TRUE,
                          ...){

 if(run_checks) check_predict(object, new_data, times, risk)

 x_new <- as.matrix(
  one_hot(x_data = new_data,
          fi = get_fctr_info(object),
          names_x_data = get_names_x(object))
 )

 if(length(times) == 1){

  return(orsf_pred_uni(object$forest, x_new, times, risk))

 }

 orsf_pred_multi(object$forest, x_new, times, risk)

}


check_new_data_names <- function(new_data,
                                 ref_names,
                                 label_new,
                                 label_ref,
                                 check_new_in_ref = FALSE,
                                 check_ref_in_new = TRUE){

 new_names <- names(new_data)

 list_new <- FALSE

 if(check_new_in_ref) list_new <- !(new_names %in% ref_names)

 list_ref <- FALSE

 if(check_ref_in_new) list_ref <- !(ref_names %in% new_names)

 error_new <- any(list_new)
 error_ref <- any(list_ref)

 if(error_new){
  out_msg_new <- paste(
   label_new, " have columns not contained in ", label_ref, ": ",
   paste_collapse(new_names[list_new], last = ' and ')
  )
 }

 if(error_ref){
  out_msg_ref <- paste(
   label_ref, " have columns not contained in ", label_new, ": ",
   paste_collapse(ref_names[list_ref], last = ' and ')
  )
 }

 if(error_new && error_ref){
  out_msg <- c(out_msg_new, '\n Also, ', out_msg_ref)
 }

 if (error_new && !error_ref) {
  out_msg <- c(out_msg_new)
 }

 if (!error_new && error_ref){
  out_msg <- c(out_msg_ref)
 }

 any_error <- error_new | error_ref

 if(any_error){
  stop(out_msg, call. = FALSE)
 }

}

check_new_data_types <- function(new_data,
                                 ref_names,
                                 ref_types,
                                 label_new,
                                 label_ref){

 var_types <- vector(mode = 'character', length = length(ref_names))

 for(i in seq_along(ref_names)){
  var_types[i] <- class(new_data[[ ref_names[i] ]])[1]
 }

 bad_types <- which(var_types != ref_types)

 if(!is_empty(bad_types)){

  vars_to_list <- ref_names[bad_types]
  types_to_list <- var_types[bad_types]

  meat <- paste0('<', vars_to_list, '> has type <',
                 types_to_list, '>', " in ", label_new,
                 "; type <", ref_types[bad_types], "> in ",
                 label_ref, collapse = '\n')

  msg <- paste("some variables in ", label_new,
               " have different type in ",
               label_ref, ":\n", meat)

  stop(msg, call. = FALSE)

 }

}

check_new_data_fctrs <- function(new_data,
                                 names_x,
                                 fi_ref,
                                 label_new){

 fctr_check(new_data, names_x)

 fi_new <- fctr_info(new_data, names_x)

 for(fi_col in fi_ref$cols){
  fctr_check_levels(ref = fi_ref$lvls[[fi_col]],
                    new = fi_new$lvls[[fi_col]],
                    name = fi_col,
                    label_ref = "training data",
                    label_new = label_new)
 }

}

fctr_check_levels <- function(ref,
                              new,
                              name,
                              label_ref,
                              label_new){

 list_new  <- !(new %in% ref)

 if(any(list_new)){

  out_msg <- paste0(
   "variable ", name, " in ", label_new,
   " has levels not contained in ", label_ref, ": ",
   paste_collapse(new[list_new], last = ' and ')
  )

  stop(out_msg, call. = FALSE)

 }


}

check_predict <- function(object, new_data, times, risk){

 if(!is.null(new_data)){

  check_arg_is(arg_value = new_data,
               arg_name = 'new_data',
               expected_class = 'data.frame')

 }

 if(!is.null(times)){

  check_arg_type(arg_value = times,
                 arg_name = 'times',
                 expected_type = 'numeric')

  check_arg_gt(arg_value = times,
               arg_name = 'times',
               bound = 0)

 }

 if(!is.null(risk)){

  check_arg_type(arg_value = risk,
                 arg_name = 'risk',
                 expected_type = 'logical')

  check_arg_length(arg_value = risk,
                   arg_name = 'risk',
                   expected_length = 1)

 }

 if(any(times > get_max_time(object))){

  stop("prediction times should ",
       "be <= max follow-up time ",
       "observed in training data: ",
       get_max_time(object),
       ". You may bypass this error by setting run_checks = FALSE",
       call. = FALSE)

 }

 if(!all(order(times) == seq(length(times)))){
  stop("times must be entered in ascending order, e.g.,",
       "times = c(5, 10) instead of times = c(10, 5)",
       call. = FALSE)
 }

 check_new_data_names(new_data  = new_data,
                      ref_names = get_names_x(object),
                      label_new = "new_data",
                      label_ref = 'training data')

 check_new_data_types(new_data  = new_data,
                      ref_names = get_names_x(object),
                      ref_types = get_types_x(object),
                      label_new = "new_data",
                      label_ref = 'training data')

 check_new_data_fctrs(new_data  = new_data,
                      names_x   = get_names_x(object),
                      fi_ref    = get_fctr_info(object),
                      label_new = "new_data")


}

