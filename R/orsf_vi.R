

#' ORSF variable importance
#'
#' @param object an object of class 'aorsf'.
#'
#' @return the `object` with variable importance attached.
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
orsf_vi <- function(object){

 if(!is_aorsf(object)) stop("object must inherit from 'aorsf' class.",
                            call. = FALSE)

 if(!has_oobag(object)){
  stop("cannot compute variable importance if the aorsf object does",
       " not have out-of-bag error (see oobag_pred in ?orsf).",
       call. = FALSE)
 }

 cstat <- last_value(object$eval_oobag$c_harrell[, 1, drop=TRUE])

 y <- object$data_train[, get_names_y(object)]

 sorted <- order(y[, 1], -y[, 2])

 x <- as.matrix(
  one_hot(x_data = object$data_train,
          fi = get_fctr_info(object),
          names_x_data = get_names_x(object))
 )

 out <- orsf_oob_vi(x = x[sorted, ],
                    cstat = cstat,
                    forest = object$forest,
                    time_pred_ = object$time_pred)

 rownames(out) <- colnames(x)

 rev(out[order(out), , drop=TRUE])

}




