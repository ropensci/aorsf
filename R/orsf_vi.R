

#' ORSF variable importance
#'
#' @param object an object of class 'aorsf'.
#'
#' @param group_factors (_logical_) if `TRUE`, the importance of factor
#'   variables will be reported overall by aggregating the importance
#'   of individual levels of the factor. If `FALSE`, the importance of
#'   individual factor levels will be returned.
#'
#' @return the `object` with variable importance attached.
#'
#' @details
#'
#' The method used to compute variable importance with ORSF is called
#'  'negation importance'. Each variable is assessed separately by
#'  multiplying the variable's coefficients by -1 and then determining
#'  how much the model's performance changes. The worse the model's
#'  performance after negating all coefficients for a given variable,
#'  the more important the variable.
#'
#' @export
#'
#' @examples
#'
#' fit <- orsf(pbc_orsf,
#'             Surv(time, status) ~ . - id,
#'             oobag_pred = TRUE)
#'
#' orsf_vi(fit)
#'
orsf_vi <- function(object, group_factors = TRUE){

 if(!is_aorsf(object)) stop("object must inherit from 'aorsf' class.",
                            call. = FALSE)

 if(!contains_oobag(object)){
  stop("cannot compute variable importance if the aorsf object does",
       " not have out-of-bag error (see oobag_pred in ?orsf).",
       call. = FALSE)
 }

 if(contains_vi(object)) return(object$importance)

 cstat <- last_value(object$eval_oobag$c_harrell[, 1, drop=TRUE])

 y <- as.matrix(object$data_train[, get_names_y(object)])

 sorted <- order(y[, 1], -y[, 2])

 x <- as.matrix(
  one_hot(x_data = object$data_train,
          fi = get_fctr_info(object),
          names_x_data = get_names_x(object))
 )

 out <- orsf_oob_vi(x = x[sorted, ],
                    y = y[sorted, ],
                    cstat = cstat,
                    forest = object$forest,
                    time_pred_ = object$time_pred)

 rownames(out) <- colnames(x)

 if(group_factors) {

  fi <- get_fctr_info(object)

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

 rev(out[order(out), , drop=TRUE])

}




