

penalized_cph <- function(x_node,
                          y_node,
                          w_node,
                          alpha,
                          df_target){

 suppressWarnings(
  fit <- try(
   glmnet::glmnet(x = x_node,
                  y = y_node,
                  weights = w_node,
                  alpha = alpha,
                  family = "cox"),
   silent = TRUE
  )
 )

 if(class(fit)[1] == 'try-error'){
  return(matrix(0, nrow=ncol(x_node), ncol=1))
 }

 for(i in seq_along(fit$df)){
  if(fit$df[i] >= df_target || i == length(fit$df)){
   return(matrix(fit$beta[, i, drop=TRUE], ncol = 1))
  }
 }

}

