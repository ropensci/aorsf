

penalized_cph <- function(x_node,
                          y_node,
                          alpha=1/2){

 # fit <- try(
 #  glmnet::glmnet(x = x_node,
 #                 y = y_node,
 #                 family = "cox",
 #                 alpha = alpha),
 #  silent = TRUE
 # )
 #
 #
 # if(class(fit)[1] == 'try-error'){
 #  return(matrix(0, nrow=ncol(x_node), ncol=1))
 # }

 fit <-
  glmnet::glmnet(x = x_node,
                 y = y_node,
                 family = "cox",
                 alpha = alpha)

 matrix(fit$beta[, ncol(fit$beta), drop=TRUE], ncol = 1)

}

