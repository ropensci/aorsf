



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
ostree_pred_leaf <- function(tree, x_new){

 parts <- rep(0, nrow(x_new))
 N <- tree$nodes

 for(nn in seq_along(N)){

  parts_in_node <- parts == names(N)[nn]

  if(any(parts_in_node)){

   i <- which(parts_in_node)

   parts[i] <- ifelse(
    x_new[i, N[[nn]]$vars] %*% N[[nn]]$vals - N[[nn]]$cutpoint < 0,
    yes = N[[nn]]$child_left,
    no = N[[nn]]$child_right
   )

  }

 }

 as.character(parts)

}

ostree_predict <- function(tree, x_new, times){



}


# x_leaves <- data.table(node = ostree_pred_leaf(tree, x_new = x_new))
# x_leaves[, node := as.character(node)]
# setkey(x_leaves, node)
#
# times <- c(50, 100, 150, 5000)
#
# dt_times <- as.data.table(
#  expand.grid(node = names(tree$leaves),
#              time = times,
#              stringsAsFactors = FALSE)
# )
#
# setkey(dt_times, node, time)
#
# tmp <- dt_leaves[dt_times,  roll = -Inf]
# tmp[, surv := nafill(surv, 'locf'), by = node]
# dcast(tmp[x_leaves, ], node ~ time)
#
#
# x_leaves











