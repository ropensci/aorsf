

#' Title
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
orsf_interaction <- function(object){

 dt_betas <- dt_means <- list()

 xnames <- get_names_x(object, one_hot_names = TRUE)

 for(t in seq_along(object$forest)){

  tree <- object$forest[[t]]

  nodes <- which(tree$children_left > 0)

  tree_means <- tree_coefs <- matrix(NA_real_,
                                     nrow = length(nodes),
                                     ncol = length(xnames))

  colnames(tree_coefs) <- colnames(tree_means) <- xnames

  for(i in seq(nrow(tree_coefs))){

   ii <- nodes[i]

   cols <- tree$col_indices[, ii] + 1

   beta_nonzero <- which(tree$betas[, ii] != 0)

   tree_coefs[i, cols[beta_nonzero]] <- tree$betas[beta_nonzero, ii]
   tree_means[i, cols[beta_nonzero]] <- tree$x_mean[beta_nonzero, ii]

  }

  dt_betas[[t]] <- as.data.table(tree_coefs)
  dt_means[[t]] <- as.data.table(tree_means)

 }

 dt_betas <- as.matrix(rbindlist(dt_betas))
 dt_means <- as.matrix(rbindlist(dt_means))

 imat <- matrix(NA_real_,
                nrow = length(xnames),
                ncol = length(xnames))

 colnames(imat) <- rownames(imat) <- xnames

 n_vars <- length(xnames)

 i_vals <- seq(n_vars)

 for(i in i_vals){

  imat[xnames[-i], xnames[i]] <-
   t(cor(dt_betas[,i], dt_means[,-i], use = 'pairwise.complete.obs')^2)

 }

 for(i in i_vals){

  for(j in seq(i, n_vars)){
   if(j > i){
    imat[i, j] <- (imat[i, j] + imat[j, i]) / 2
    imat[j, i] <- NA_real_
   }

  }

 }

 dt <- as.data.table(imat, keep.rownames = 'v1')
 dt <- melt(dt, id.vars = 'v1', variable.name = 'v2', na.rm = TRUE)
 dt <- dt[order(value, decreasing = TRUE), ]
 dt[, perc_drop := shift(value, n=1, fill=0) / value]
 dt[]

 dt

}
