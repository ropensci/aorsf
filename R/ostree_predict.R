



#' Title
#'
#' @param tree
#' @param x_new
#'
#' @return
#' @export
#'
#' @examples
#'
#' x_new <- flchain_x[1:10, ]
#'


ostree_predict <- function(tree, x_new, times){

 leaf_nodes = ostree_pred_leaf(x_new,
                               tree$betas,
                               tree$col_indices,
                               tree$cut_points,
                               tree$children_left)

 1-ostree_pred_surv(x_new, tree$leaves, leaf_nodes, times)


}

# TODO: throw an error if any(times > max(tree$times))
