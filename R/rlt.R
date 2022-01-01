




rlt <- function(x_node,
                y_node,
                w_node,
                n_keep,
                n_tree){


 rlt_vars <- seq(ncol(x_node))

 out <- matrix(0, nrow = ncol(x_node), ncol = 1)

 if(n_keep < nrow(out)){

  colnames(x_node) <- paste0('x', seq(ncol(x_node)))

  suppressWarnings(
   fit <- try(
    randomForestSRC::rfsrc.fast(
     Surv(time, status) ~ .,
     data = as.data.frame(cbind(y_node, x_node)),
     ntree = n_tree,
     nodesize = 3,
     importance = 'anti'
    ),
    silent = FALSE
   )
  )

  if(class(fit)[1] != 'try-error'){

   if(any(fit$importance > 0)) {

    sorted_index <- order(fit$importance, decreasing = TRUE)

    rlt_vars <- sorted_index[seq(n_keep)]

   }

  }

 }

 rlt_coefs <- penalized_cph(x_node = x_node[, rlt_vars],
                            y_node = y_node,
                            w_node = w_node,
                            alpha = 1/2,
                            df_target = length(rlt_vars))

 if(any(rlt_coefs != 0)){

  out[rlt_vars, ] <- rlt_coefs

 }

 out

}
