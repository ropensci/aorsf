

#' Variable Importance
#'
#' Estimate the importance of individual predictor variables using
#'  oblique random forests.
#'
#' @inheritParams predict.ObliqueForest
#'
#' @param group_factors (_logical_) `r roxy_group_factors()`
#'
#' @param importance `r roxy_importance_header()`
#' - `r roxy_importance_anova()`
#' - `r roxy_importance_negate()`
#' - `r roxy_importance_permute()`
#'
#' @param oobag_fun `r roxy_oobag_fun_header()` after negating coefficients
#'   (if importance = 'negate') or permuting the values of a predictor
#'   (if importance = 'permute')
#' - `r roxy_oobag_fun_default()`
#' - `r roxy_oobag_fun_user()`
#'     + `r roxy_oobag_fun_inputs()`
#'     + `r roxy_oobag_fun_ymat()`
#'     + `r roxy_oobag_fun_svec()`
#'     + `r roxy_oobag_fun_return()`
#'     + the same `oobag_fun` should have been used when you created `object`
#'       so that the initial value of out-of-bag prediction accuracy is
#'       consistent with the values that will be computed while variable
#'       importance is estimated.
#'
#' For more details, see the out-of-bag
#' [vignette](https://docs.ropensci.org/aorsf/articles/oobag.html).
#'
#' @section Variable importance methods:
#'
#' __negation importance__: `r roxy_vi_describe('negate')`
#'
#' __permutation importance__: `r roxy_vi_describe('permute')`
#'
#' __analysis of variance (ANOVA) importance__: `r roxy_vi_describe('anova')`
#'
#' @return `orsf_vi` functions return a named numeric vector.
#'
#' - Names of the vector are the predictor variables used by `object`
#' - Values of the vector are the estimated importance of the given predictor.
#'
#' The returned vector is sorted from highest to lowest value, with higher
#'  values indicating higher importance.
#'
#' @details
#'
#' When an `ObliqueForest` object is grown with importance = 'anova',
#'  'negate', or 'permute', the output will have a vector of importance
#'  values based on the requested type of importance. However, `orsf_vi()`
#'  can be used to compute variable importance after growing a forest
#'  or to compute a different type of importance.
#'
#' `orsf_vi()` is a general purpose function to extract or compute variable
#'  importance estimates from an `ObliqueForest` object (see [orsf]).
#'  `orsf_vi_negate()`, `orsf_vi_permute()`, and `orsf_vi_anova()` are wrappers
#'  for `orsf_vi()`. The way these functions work depends on whether the
#'  `object` they are given already has variable importance estimates in it
#'  or not (see examples).
#'
#'
#' @export
#'
#' @includeRmd Rmd/orsf_vi_examples.Rmd
#'
#' @references
#'
#' 1. `r cite("harrell_1982")`
#' 1. `r cite("breiman_2001")`
#' 1. `r cite("menze_2011")`
#' 1. `r cite("jaeger_2022")`
#'
#'
orsf_vi <- function(object,
                    group_factors = TRUE,
                    importance = NULL,
                    oobag_fun = NULL,
                    n_thread = NULL,
                    verbose_progress = NULL,
                    ...){

 check_dots(list(...), .f = orsf_vi)

 orsf_vi_(object,
          group_factors = group_factors,
          importance = importance,
          oobag_fun = oobag_fun,
          n_thread = n_thread,
          verbose_progress = verbose_progress)


}

#' @rdname orsf_vi
#' @export
orsf_vi_negate <-
 function(object,
          group_factors = TRUE,
          oobag_fun = NULL,
          n_thread = NULL,
          verbose_progress = NULL,
          ...) {
  check_dots(list(...), .f = orsf_vi_negate)
  orsf_vi_(object,
           group_factors,
           importance = 'negate',
           oobag_fun = oobag_fun,
           n_thread = n_thread,
           verbose_progress = verbose_progress)
 }

#' @rdname orsf_vi
#' @export
orsf_vi_permute <-
 function(object,
          group_factors = TRUE,
          oobag_fun = NULL,
          n_thread = NULL,
          verbose_progress = NULL,
          ...) {
  check_dots(list(...), .f = orsf_vi_permute)
  orsf_vi_(object,
           group_factors,
           importance = 'permute',
           oobag_fun = oobag_fun,
           n_thread = n_thread,
           verbose_progress = verbose_progress)
 }

#' @rdname orsf_vi
#' @export
orsf_vi_anova <- function(object,
                          group_factors = TRUE,
                          verbose_progress = NULL,
                          ...) {
 check_dots(list(...), .f = orsf_vi_anova)
 orsf_vi_(object,
          group_factors,
          importance = 'anova',
          oobag_fun = NULL,
          n_thread = 0,
          verbose_progress = FALSE)
}

#' Variable importance working function
#'
#' @inheritParams orsf_vi_negate
#' @param importance the type of variable importance technique to use.
#'
#' @noRd
#'
orsf_vi_ <- function(object,
                     group_factors,
                     importance,
                     oobag_fun,
                     n_thread,
                     verbose_progress){

 check_arg_is(object, arg_name = 'object', expected_class = 'ObliqueForest')

 if(is.null(object$data)){
  stop("training data were not found in object, ",
       "but are needed to collapse factor importance values. ",
       "Did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)
 }

 type_vi <- object$importance_type

 if(type_vi == 'none' && is.null(importance))
  stop("object has no variable importance values to extract and the type ",
       "of importance to compute is not specified. ",
       "Try setting importance = 'permute', or 'negate', or 'anova' ",
       "in orsf() or setting importance = 'permute' or 'negate' in ",
       "orsf_vi().",
       call. = FALSE)

 # check n_thread

 if(!is.null(n_thread)){

  check_arg_type(arg_value = n_thread,
                 arg_name = 'n_thread',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = n_thread, arg_name = 'n_thread')

  check_arg_gteq(arg_value = n_thread, arg_name = 'n_thread', bound = 0)

  check_arg_length(arg_value = n_thread,
                   arg_name = 'n_thread',
                   expected_length = 1)

 }

 # check verbose progress

 if(!is.null(verbose_progress)){

  check_arg_type(arg_value = verbose_progress,
                 arg_name = 'verbose_progress',
                 expected_type = 'logical')

  check_arg_length(arg_value = verbose_progress,
                   arg_name = 'verbose_progress',
                   expected_length = 1)

 }


 if(!is.null(importance)){

  # check importance
  check_arg_type(arg_value = importance,
                 arg_name = 'importance',
                 expected_type = 'character')

  check_arg_length(arg_value = importance,
                   arg_name = 'importance',
                   expected_length = 1)

  check_arg_is_valid(arg_value = importance,
                     arg_name = 'importance',
                     valid_options = c("none",
                                       "anova",
                                       "negate",
                                       "permute"))

  # would someone call orsf_vi(importance = 'none')? just in case...
  if(importance == 'none'){
   warning("orsf_vi was called with importance = 'none'. Returning NULL")
   return(NULL)
  }

  if(object$importance_type != 'anova' && importance == 'anova')
   stop("ANOVA importance can only be computed while an orsf object",
        " is being fitted. To get ANOVA importance values, train your",
        " orsf object with importance = 'anova'",
        call. = FALSE)

  type_vi <- importance

 }

 if(object$sample_fraction > 0.90 && type_vi %in% c("negate", "permute")){
  stop("sample_fraction must be <= 9/10 for out-of-bag importance estimates",
       call. = FALSE)
 }

 out <- switch(
  type_vi,
  'anova' = object$get_importance_raw(),
  'negate' = orsf_vi_oobag_(object, type_vi, oobag_fun,
                            n_thread, verbose_progress),
  'permute' = orsf_vi_oobag_(object, type_vi, oobag_fun,
                             n_thread, verbose_progress)
 )

 object$get_importance_clean(out, group_factors)

}

#' Variable importance oobag working function
#'
#' used to pass data and function specification to C++
#'
#' @inheritParams orsf_vi_
#'
#' @noRd
#'
orsf_vi_oobag_ <- function(object,
                           type_vi,
                           oobag_fun,
                           n_thread,
                           verbose_progress){

 if(contains_vi(object) && is.null(oobag_fun) &&
    object$importance_type == type_vi){

  out <- object$get_importance_raw()

  return(out)

 }

 object$compute_vi(type_vi, oobag_fun, n_thread, verbose_progress)


}





