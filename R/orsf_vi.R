

#' ORSF variable importance
#'
#' Estimate the importance of individual variables using oblique random
#'   survival forests.
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
#' When an `orsf_fit` object is fitted with importance = 'anova', 'negate', or
#'  'permute', the output will have a vector of importance values based on
#'  the requested type of importance. However, you may still want to call
#'  `orsf_vi()` on this output if you want to group factor levels into one
#'  overall importance value.
#'
#' `orsf_vi()` is a general purpose function to extract or compute variable
#'   importance estimates from an `'orsf_fit'` object (see [orsf]).
#'   `orsf_vi_negate()`, `orsf_vi_permute()`, and `orsf_vi_anova()` are wrappers
#'   for `orsf_vi()`. The way these functions work depends on whether the
#'   `object` they are given already has variable importance estimates in it
#'   or not (see examples).
#'
#'
#' @export
#'
#' @includeRmd Rmd/orsf_vi_examples.Rmd
#'
#' @references
#'
#'
#' `r roxy_cite_harrell_1982()`
#'
#' `r roxy_cite_breiman_2001()`
#'
#' `r roxy_cite_menze_2011()`
#'
#' `r roxy_cite_jaeger_2023()`
#'
#'
orsf_vi <- function(object,
                    group_factors = TRUE,
                    importance = NULL,
                    oobag_fun = NULL,
                    n_thread = 1,
                    verbose_progress = FALSE,
                    ...){

 check_dots(list(...), .f = orsf_vi)

 # not sure if anyone would ever call orsf_vi(importance = 'none'),
 # but this is here just for them.
 if(!is.null(importance)){
  if(importance == 'none') importance <- NULL
 }

 if(!is.null(importance)){

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


 }

 type_vi <- object$importance_type

 if(type_vi == 'none' && is.null(importance))
  stop("object has no variable importance values to extract and the type ",
       "of importance to compute is not specified. ",
       "Try setting importance = 'permute', or 'negate', or 'anova' ",
       "in orsf() or setting importance = 'permute' or 'negate' in ",
       "orsf_vi().",
       call. = FALSE)

 if(!is.null(importance)) type_vi <- importance

 orsf_vi_(object,
          group_factors = group_factors,
          type_vi = type_vi,
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
          n_thread = 1,
          verbose_progress = FALSE,
          ...) {
  check_dots(list(...), .f = orsf_vi_negate)
  orsf_vi_(object,
           group_factors,
           type_vi = 'negate',
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
          n_thread = 1,
          verbose_progress = FALSE,
          ...) {
  check_dots(list(...), .f = orsf_vi_permute)
  orsf_vi_(object,
           group_factors,
           type_vi = 'permute',
           oobag_fun = oobag_fun,
           n_thread = n_thread,
           verbose_progress = verbose_progress)
 }

#' @rdname orsf_vi
#' @export
orsf_vi_anova <- function(object,
                          group_factors = TRUE,
                          ...) {
 check_dots(list(...), .f = orsf_vi_anova)
 orsf_vi_(object,
          group_factors,
          type_vi = 'anova',
          oobag_fun = NULL,
          verbose_progress = FALSE)
}

#' Variable importance working function
#'
#' @inheritParams orsf_vi_negate
#' @param type_vi the type of variable importance technique to use.
#'
#' @noRd
#'
orsf_vi_ <- function(object,
                     group_factors,
                     type_vi,
                     oobag_fun,
                     n_thread,
                     verbose_progress){

 if(!is_aorsf(object)) stop("object must inherit from 'ObliqueForest' class.",
                            call. = FALSE)

 if(object$importance_type != 'anova' && type_vi == 'anova')
  stop("ANOVA importance can only be computed while an orsf object",
       " is being fitted. To get ANOVA importance values, train your",
       " orsf object with importance = 'anova'",
       call. = FALSE)

 if(is.null(object$data)){
  stop("training data were not found in object, ",
       "but are needed to collapse factor importance values. ",
       "Did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)
 }

 out <- switch(
  type_vi,
  'anova' = object$get_importance_raw(),
  'negate' = orsf_vi_oobag_(object, type_vi, oobag_fun,
                            n_thread, verbose_progress),
  'permute' = orsf_vi_oobag_(object, type_vi, oobag_fun,
                             n_thread, verbose_progress)
 )

 # nan's occur if a variable was never used, so:
 out[is.nan(out)] <- 0

 if(group_factors) {

  fi <- object$get_fctr_info()

  if(!is_empty(fi$cols)){

   for(f in fi$cols[!fi$ordr]){

    f_lvls <- fi$lvls[[f]]
    f_rows <- match(paste(f, f_lvls[-1], sep = '_'), rownames(out))
    f_wts <- 1

    if(length(f_lvls) > 2){
     f_wts <- prop.table(x = table(object$data[[f]])[-1])
    }

    f_vi <- sum(out[f_rows] * f_wts, na.rm = TRUE)

    out[f_rows] <- f_vi
    rownames(out)[f_rows] <- f

   }

   if(!is_empty(fi$cols[!fi$ordr])) {
    # take extra care in case there are duplicate vi values.
    out <- data.frame(variable = rownames(out), value = out)
    out <- unique(out)
    out$variable <- NULL
    out <- as.matrix(out)
    colnames(out) <- NULL
   }



  }

 }

 rev(out[order(out), , drop=TRUE])

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





