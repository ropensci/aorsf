
#' Update Forest Parameters
#'
#' @param object `r roxy_describe_ObliqueForest(trained = FALSE)`.
#'
#' @param ... arguments to plug into [orsf] that will be used to define the
#'  update (see examples).
#'
#' @param modify_in_place (*logical*) if `TRUE`, `object` will be modified
#'   by the inputs specified in `...`. Be cautious, as modification in place
#'   will overwrite existing data. If `FALSE` (the default), `object` will
#'   be copied and then the modifications will be applied to the copy,
#'   leaving the original `object` unmodified.
#'
#' @inheritParams orsf
#'
#' @details
#'
#'  There are several dynamic inputs in `orsf` with default values of `NULL`.
#'  Specifically, these inputs are `control`, `weights`, `mtry`, `split_rule`,
#'  `split_min_stat`, `pred_type`, `pred_horizon`, `oobag_eval_function`,
#'  `tree_seeds`, and `oobag_eval_every`. If no explicit value is given for
#'  these inputs in the call, they *will be re-formed*. For example, if
#'  an initial forest includes 17 predictors, the default `mtry` is the
#'  smallest integer that is greater than or equal to the square root of 17,
#'  i.e., 5. Then, if you make an updated forest with 1 less predictor and
#'  you do not explicitly say `mtry = 5`, then `mtry` will be re-initialized
#'  in the update based on the available 16 predictors, and the resulting
#'  value of `mtry` will be 4. This is done to avoid many potential errors
#'  that would occur if the dynamic outputs were not re-initialized.
#'
#' @return an `ObliqueForest` object.
#'
#' @export
#'
#' @examples
#'
#' # initial fit has mtry of 5
#' fit <- orsf(pbc_orsf, time + status ~ . -id)
#'
#' # note that mtry is now 4 (see details)
#' fit_new <- orsf_update(fit, formula = . ~ . - edema, n_tree = 100)
#'
#' # prevent dynamic updates by specifying inputs you want to freeze.
#' fit_newer <- orsf_update(fit_new, mtry = 2)
#'
#'
orsf_update <- function(object,
                        ...,
                        modify_in_place = FALSE,
                        no_fit = NULL){

 check_arg_is(object, arg_name = 'object', expected_class = 'ObliqueForest')

 .dots <- list(...)
 user_args <- names(.dots)
 orsf_args <- setdiff(names(formals(orsf)), "...")
 user_args_unmatched <- setdiff(user_args, orsf_args)
 check_dots(.dots[user_args_unmatched], .f = orsf)

 no_fit <- no_fit %||% !object$trained

 if(modify_in_place){

  object_new <- object

 } else {

  object_new <- object$clone(deep = TRUE)

 }

 object_new$update(...)

 if(no_fit){

  object_new$untrain()

 } else {

  object_new$train()

 }

 if(modify_in_place) return(invisible(object_new))

 object_new

}





