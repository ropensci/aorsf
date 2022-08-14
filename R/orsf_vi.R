

#' ORSF variable importance
#'
#' Estimate the importance of individual variables using oblique random
#'   survival forests.
#'
#' @inheritParams predict.orsf_fit
#'
#' @param group_factors (_logical_) if `TRUE`, the importance of factor
#'   variables will be reported overall by aggregating the importance
#'   of individual levels of the factor. If `FALSE`, the importance of
#'   individual factor levels will be returned.
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
#' [vignette](https://bcjaeger.github.io/aorsf/articles/oobag.html).
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
#' @examples
#'
#' # first example ----------------------------------------------------------
#'
#' # fit an ORSF using default values, and get the default vi (anova)
#'
#' fit_default <- orsf(pbc_orsf,
#'                     Surv(time, status) ~ . - id)
#'
#' # the printed output will indicate the type of vi used
#'
#' fit_default
#'
#' # the 'raw' vi values are stored in fit_default:
#'
#' fit_default$importance
#'
#' # these are 'raw' because the vi values for factors have not been
#' # aggregated into a single vi value. Currently there is one vi value
#' # for k-1 levels of a k level factor. For example, you can see edema_1
#' # and edema_0.5 in the importance values above because edema is a factor
#' # variable with levels of 0, 0.5, and 1. To get aggregated values of vi
#' # across all levels of each factor, just call orsf_vi with group_factors
#' # set to TRUE (the default)
#'
#' orsf_vi(fit_default)
#'
#' # orsf_vi knows that fit_default was fit using anova vi
#' # to verify this, see that orsf_vi and orsf_vi_anova are the same
#'
#' orsf_vi_anova(fit_default)
#'
#'
#' # second example ---------------------------------------------------------
#'
#' # fit an ORSF without vi, then add vi later
#'
#' fit_no_vi <- orsf(pbc_orsf,
#'                   Surv(time, status) ~ . - id,
#'                   importance = 'none')
#'
#' # Note: you can't call orsf_vi_anova() on fit_no_vi because anova
#' # VI can only be computed while the forest is being grown.
#'
#' orsf_vi_negate(fit_no_vi)
#'
#' orsf_vi_permute(fit_no_vi)
#'
#' # third example ----------------------------------------------------------
#'
#' # fit an ORSF and compute vi at the same time
#'
#' fit_permute_vi <- orsf(pbc_orsf,
#'                        Surv(time, status) ~ . - id,
#'                        importance = 'permute')
#'
#' # get the vi instantly (i.e., it doesn't need to be computed again)
#'
#' orsf_vi_permute(fit_permute_vi)
#'
#' # You can still get negation vi from this fit, but it needs to be computed
#'
#' orsf_vi_negate(fit_permute_vi)
#'
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
#' `r roxy_cite_jaeger_2022()`
#'
#'
orsf_vi <- function(object,
                    group_factors = TRUE,
                    importance = NULL,
                    oobag_fun = NULL,
                    ...){

 check_dots(list(...), .f = orsf_vi)

 # not sure if anyone would ever call orsf_vi(importance = 'none'),
 # but this is here just for them.
 if(!is.null(importance)){
  if(importance == 'none') importance <- NULL
 }

 check_orsf_inputs(importance = importance)

 type_vi <- get_importance(object)

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
          oobag_fun = oobag_fun)


}

#' @rdname orsf_vi
#' @export
orsf_vi_negate <- function(object, group_factors = TRUE, oobag_fun = NULL, ...){
 check_dots(list(...), .f = orsf_vi_negate)
 orsf_vi_(object, group_factors, type_vi = 'negate', oobag_fun = oobag_fun)
}

#' @rdname orsf_vi
#' @export
orsf_vi_permute <- function(object, group_factors = TRUE, oobag_fun = NULL, ...){
 check_dots(list(...), .f = orsf_vi_permute)
 orsf_vi_(object, group_factors, type_vi = 'permute', oobag_fun = oobag_fun)
}

#' @rdname orsf_vi
#' @export
orsf_vi_anova <- function(object, group_factors = TRUE, ...){
 check_dots(list(...), .f = orsf_vi_anova)
 orsf_vi_(object, group_factors, type_vi = 'anova', oobag_fun = NULL)
}

#' Variable importance working function
#'
#' @inheritParams orsf_vi_negate
#' @param type_vi the type of variable importance technique to use.
#'
#' @noRd
#'
orsf_vi_ <- function(object, group_factors, type_vi, oobag_fun = NULL){

 #' @srrstats {G2.8} *As part of initial pre-processing, run checks on inputs to ensure that all other sub-functions receive inputs of a single defined class or type.*

 if(!is_aorsf(object)) stop("object must inherit from 'orsf_fit' class.",
                            call. = FALSE)

 if(get_importance(object) != 'anova' && type_vi == 'anova')
  stop("ANOVA importance can only be computed while an orsf object",
       " is being fitted. To get ANOVA importance values, re-grow your",
       " orsf object with importance = 'anova'",
       call. = FALSE)

 out <- switch(type_vi,
               'anova' = as.matrix(object$importance),
               'negate' = orsf_vi_oobag_(object, type_vi, oobag_fun),
               'permute' = orsf_vi_oobag_(object, type_vi, oobag_fun))

 if(group_factors) {

  if(is.null(object$data)){
   stop("training data were not found in object, ",
        "but are needed to collapse factor importance values. ",
        "Did you use attach_data = FALSE when ",
        "running orsf()?", call. = FALSE)
  }

  fi <- get_fctr_info(object)

  if(!is_empty(fi$cols)){

   for(f in fi$cols[!fi$ordr]){

    f_lvls <- fi$lvls[[f]]
    f_rows <- match(paste(f, f_lvls[-1], sep = '_'), rownames(out))
    f_wts <- 1

    if(length(f_lvls) > 2){
     f_wts <- prop.table(x = table(object$data[[f]])[-1])
    }

    f_vi <- sum(out[f_rows] * f_wts)

    out[f_rows] <- f_vi
    rownames(out)[f_rows] <- f

   }

   if(!is_empty(fi$cols[!fi$ordr])) out <- unique(out)

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
orsf_vi_oobag_ <- function(object, type_vi, oobag_fun){

 if(!contains_oobag(object)){
  stop("cannot compute ",
       switch(type_vi, 'negate' = 'negation', 'permute' = 'permutation'),
       " importance if the orsf_fit object does not have out-of-bag error",
       " (see oobag_pred in ?orsf).",
       call. = FALSE)
 }

 if(contains_vi(object) &&
    is.null(oobag_fun) &&
    get_importance(object) == type_vi){

  out <- matrix(object$importance, ncol = 1)

  rownames(out) <- names(object$importance)

  return(out)

 }

 if(is.null(oobag_fun)){

  f_oobag_eval <- function(x) x
  type_oobag_eval <- 'H'

 } else {

  check_oobag_fun(oobag_fun)
  f_oobag_eval <- oobag_fun
  type_oobag_eval <- 'U'

 }

 y <- as.matrix(object$data[, get_names_y(object)])

 # Put data in the same order that it was in when object was fit
 sorted <- order(y[, 1], -y[, 2])

 x <- as.matrix(
  ref_code(x_data = object$data,
           fi = get_fctr_info(object),
           names_x_data = get_names_x(object))
 )

 if(is.null(oobag_fun)) {

  last_eval_stat <-
   last_value(object$eval_oobag$stat_values[, 1, drop=TRUE])

 } else {

  last_eval_stat <-
   f_oobag_eval(y_mat = y, s_vec = object$surv_oobag)

 }

 f_oobag_vi <- switch(
  type_vi,
  'negate' = orsf_oob_negate_vi,
  'permute' = orsf_oob_permute_vi
 )

 out <- f_oobag_vi(x = x[sorted, ],
                   y = y[sorted, ],
                   forest = object$forest,
                   last_eval_stat = last_eval_stat,
                   time_pred_ = object$pred_horizon,
                   f_oobag_eval = f_oobag_eval,
                   type_oobag_eval_ = type_oobag_eval)

 rownames(out) <- colnames(x)

 out

}





