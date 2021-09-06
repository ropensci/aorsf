



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



}

# TODO: throw an error if any(times > max(tree$times))


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











