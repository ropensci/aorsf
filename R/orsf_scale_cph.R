


#' Scale input data
#'
#' These functions are exported so that users may access internal routines
#'   that are used to scale inputs when [orsf_control_cph] is used.
#'
#'
#' @param x_mat (_numeric matrix_) a matrix with values to be scaled or
#'   unscaled. Note that `orsf_unscale_cph` will only accept `x_mat` inputs
#'   that have an attribute containing transform values, which are added
#'   automatically by `orsf_scale_cph`.
#'
#' @param w_vec (_numeric vector_) an optional vector of weights. If no weights
#'   are supplied (the default), all observations will be equally weighted. If
#'   supplied, `w_vec` must have length equal to `nrow(x_mat)`.
#'
#' @return the scaled or unscaled `x_mat`.
#'
#' @details The data are transformed by first subtracting the mean and then
#'   multiplying by the scale. An inverse transform can be completed using
#'   `orsf_unscale_cph` or by dividing each column by the corresponding scale
#'   and then adding the mean.
#'
#'   The values of means and scales are stored in an attribute of the output
#'   returned by `orsf_scale_cph` (see examples)
#'
#' @export
#'
#' @examples
#'
#' x_mat <- as.matrix(pbc_orsf[, c('bili', 'age', 'protime')])
#'
#' head(x_mat)
#'
#' x_scaled <- orsf_scale_cph(x_mat)
#'
#' head(x_scaled)
#'
#' attributes(x_scaled) # note the transforms attribute
#'
#' x_unscaled <- orsf_unscale_cph(x_scaled)
#'
#' head(x_unscaled)
#'
#' # numeric difference in x_mat and x_unscaled should be practically 0
#' max(abs(x_mat - x_unscaled))

orsf_scale_cph <- function(x_mat, w_vec = NULL){

 check_arg_is(arg_value = x_mat,
              arg_name = 'x_mat',
              expected_class = 'matrix')

 check_arg_type(arg_value = x_mat,
                arg_name = 'x_mat',
                expected_type = 'numeric')

 if(is_empty(x_mat))
  stop("x_mat is empty", call. = FALSE)

 if(is.null(w_vec))
  w_vec <- rep(1, nrow(x_mat))

 check_arg_type(arg_value = w_vec,
                arg_name = 'w_vec',
                expected_type = 'numeric')

 check_arg_gt(arg_value = w_vec,
              arg_name = 'w_vec',
              bound = 0)

 if(length(w_vec) != nrow(x_mat))
  stop("w_vec must have length equal to the number of rows in x_mat",
       call. = FALSE)

 # pass x[, ] instead of x to prevent x from being modified in place.
 output <- cph_scale(x_mat[, ], w_vec)

 colnames(output$x_scaled) <- colnames(x_mat)
 colnames(output$x_transforms) <- c("mean", "scale")

 out <- output$x_scaled
 attr(out, 'transforms') <- output$x_transforms

 out

}

#' @rdname orsf_scale_cph
#' @export
orsf_unscale_cph <- function(x_mat){

 if(is.null(attr(x_mat, 'transforms')))
  stop('x_mat does not have the \'transforms\' attribute, ',
       'which is needed to unscale x_mat',
       call. = FALSE)

 transforms <- attr(x_mat, 'transforms')

 # unnecessary but conceptually helpful assignment
 out <- x_mat

 for(i in seq(ncol(out))){

  out[, i] <- out[, i] / transforms[i, 'scale'] + transforms[i, 'mean']

 }

 out

}



