

#' Title
#'
#' @param object
#' @param x_new
#'
#' @return
#' @export
#'
#' @examples
#'
#' x_new <- flchain_x[1:10, ]
#'

# TODO: throw an error if any(times > max(tree$times))

orsf_predict <- function(object, x_new, times){


 x_new_scale_cph(x_new, object$x_transforms)



 for(tree in seq_along(object$forest)){

  ostree_predict_(object$forest[[tree]]$betas,
                  object$forest[[tree]]$col_indices,
                  object$forest[[tree]]$cut_points,
                  object$forest[[tree]]$children_left,
                  object$forest[[tree]]$leaf_nodes,
                  x_new,
                  times)

 }

 x_new_unscale_cph(x_new, object$x_transforms)


}
