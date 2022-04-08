

#' ORSF variable importance
#'
#' Determine the importance of individual variables using 'negation
#'   importance.' See 'Details' for definition of negation importance.
#'
#' @param object an object of class 'aorsf'.
#'
#' @param group_factors (_logical_) if `TRUE`, the importance of factor
#'   variables will be reported overall by aggregating the importance
#'   of individual levels of the factor. If `FALSE`, the importance of
#'   individual factor levels will be returned.
#'
#' @inheritParams orsf
#'
#' @return a named vector. Names indicate predictors, values indicate importance.
#'  The vector is sorted from highest to lowest value, with higher values
#'  indicating higher importance.
#'
#' @details
#'
#' __negation importance__: Each variable is assessed separately by
#'  multiplying the variable's coefficients by -1 and then determining
#'  how much the model's performance changes. The worse the model's
#'  performance after negating coefficients for a given variable,
#'  the more important the variable.
#'
#' __ANOVA importance__: ANOVA importance computes a p-value for each
#' coefficient in each linear combination of variables in each decision
#' tree of an oRF. Following the definition proposed by Menze et al,
#' ANOVA importance in aorsf for an individual variable is the proportion
#' of times a p-value for its coefficient is < 0.10.
#'
#' __Disclaimer__: Negation importance is currently in development and its routine may be tweaked in future updates. ANOVA importance has been published by Menze et al. and is in a more stable lifecycle than negation importance.
#'
#' @export
#'
#' @examples
#'
#' fit <- orsf(pbc_orsf,
#'             Surv(time, status) ~ . - id,
#'             oobag_pred = TRUE)
#'
#' orsf_vi_negate(fit)
#' orsf_vi_anova(fit)
#'
#' @references
#'
#' Menze, Bjoern H., et al. On oblique random forests.
#' *Joint European Conference on Machine Learning and Knowledge Discovery in Databases*.
#'  Springer, Berlin, Heidelberg, 2011.
#'  DOI: 10.1007/978-3-642-23783-6_29
#'
orsf_vi_negate <- function(object, group_factors = TRUE, oobag_fun = NULL){
 orsf_vi_(object, group_factors, type_vi = 'negate', oobag_fun = oobag_fun)
}


#' @rdname orsf_vi_negate
#' @export
orsf_vi_anova <- function(object, group_factors = TRUE){
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

 if(!is_aorsf(object)) stop("object must inherit from 'aorsf' class.",
                            call. = FALSE)


 switch(type_vi,

  'anova' = {
   out <- object$signif_means
  },

  'negate' = {

   if(!contains_oobag(object)){
    stop("cannot compute negation importance if the aorsf object does",
         " not have out-of-bag error (see oobag_pred in ?orsf).",
         call. = FALSE)
   }

   if(contains_vi(object) && is.null(oobag_fun)){

    out <- matrix(object$importance, ncol = 1)
    rownames(out) <- names(object$importance)


   } else {

    if(is.null(oobag_fun)){

     f_oobag_eval <- function(x) x
     type_oobag_eval <- 'H'

    } else {

     check_oobag_fun(oobag_fun)
     f_oobag_eval <- oobag_fun
     type_oobag_eval <- 'U'

    }

    y <- as.matrix(object$data_train[, get_names_y(object)])

    # Put data in the same order that it was in when object was fit
    sorted <- order(y[, 1], -y[, 2])

    x <- as.matrix(
     ref_code(x_data = object$data_train,
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


    out <- orsf_oob_vi(x = x[sorted, ],
                       y = y[sorted, ],
                       last_eval_stat = last_eval_stat,
                       forest = object$forest,
                       time_pred_ = object$pred_horizon,
                       f_oobag_eval = f_oobag_eval,
                       type_oobag_eval_ = type_oobag_eval)

    rownames(out) <- colnames(x)

   }

  }

 )

 if(group_factors) {

  fi <- get_fctr_info(object)

  if(!is_empty(fi$cols)){

   for(f in fi$cols[!fi$ordr]){

    f_lvls <- fi$lvls[[f]]
    f_rows <- match(paste(f, f_lvls[-1], sep = '_'), rownames(out))
    f_wts <- 1

    if(length(f_lvls) > 2){
     # browser()
     f_wts <- prop.table(x = table(object$data_train[[f]])[-1])
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







