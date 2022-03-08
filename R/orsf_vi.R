

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
orsf_vi_negate <- function(object, group_factors = TRUE){
 orsf_vi_(object, group_factors, type = 'negate')
}


#' @rdname orsf_vi_negate
#' @export
orsf_vi_anova <- function(object, group_factors = TRUE){
 orsf_vi_(object, group_factors, type = 'anova')
}

#' Variable importance working function
#'
#' @inheritParams orsf_vi_negate
#' @param type the type of variable selection technique to use.
#'
#' @noRd
#'
orsf_vi_ <- function(object, group_factors, type){

 #' @srrstats {G2.8} *As part of initial pre-processing, run checks on inputs to ensure that all other sub-functions receive inputs of a single defined class or type.*

 if(!is_aorsf(object)) stop("object must inherit from 'aorsf' class.",
                            call. = FALSE)


 switch(type,

  'anova' = {
   out <- object$signif_means
   rownames(out) <- get_names_x(object, ref_code_names = TRUE)
  },

  'negate' = {

   if(!contains_oobag(object)){
    stop("cannot compute negation importance if the aorsf object does",
         " not have out-of-bag error (see oobag_pred in ?orsf).",
         call. = FALSE)
   }

   if(contains_vi(object)){

    out <- object$importance

   } else {

    cstat <- last_value(object$eval_oobag$c_harrell[, 1, drop=TRUE])

    y <- as.matrix(object$data_train[, get_names_y(object)])

    sorted <- order(y[, 1], -y[, 2])

    x <- as.matrix(
     ref_code(x_data = object$data_train,
              fi = get_fctr_info(object),
              names_x_data = get_names_x(object))
    )

    out <- orsf_oob_vi(x = x[sorted, ],
                       y = y[sorted, ],
                       cstat = cstat,
                       forest = object$forest,
                       time_pred_ = object$pred_horizon)

    rownames(out) <- colnames(x)

   }

  }

 )

 if(group_factors) {

  fi <- get_fctr_info(object)

  if(!is_empty(fi$cols)){

   for(f in fi$cols[!fi$ordr]){

    f_lvls <- fi$lvls[[f]]
    f_rows <- which(rownames(out) %in% paste(f, f_lvls[-1], sep = '_'))
    f_wts <- 1

    if(length(f_lvls) > 2)
     f_wts <- prop.table(x = table(object$data_train[[f]])[-1])

    f_vi <- sum(out[f_rows] * f_wts)

    out[f_rows] <- f_vi
    rownames(out)[f_rows] <- f

   }

   if(!is_empty(fi$cols[!fi$ordr])) out <- unique(out)

  }

 }

 rev(out[order(out), , drop=TRUE])

}



orsf_vi_menze <- function(object){

 vi <- as.numeric(object$signif_means)
 names(vi) <- get_names_x(object, ref_code_names = TRUE)
 sort(vi)

}





