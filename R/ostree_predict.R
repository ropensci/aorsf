



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

# TODO: throw an error if any(times > max(object$times))

ostree_predict <- function(object, x_new, x_transforms, times){

 ostree_predict_(object$betas,
                 object$col_indices,
                 object$cut_points,
                 object$children_left,
                 object$leaf_nodes,
                 x_new[, ], # prevent modification in place
                 x_transforms,
                 times)

}


