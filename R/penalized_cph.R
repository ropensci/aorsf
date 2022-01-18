

#' Call glmnet from orsf()
#'
#' @param x_node x-matrix for the current node
#' @param y_node y-matrix for the current node
#' @param w_node weights vector for the current node
#' @param alpha see `alpha` in [glmnet::glmnet]
#' @param df_target see `dfmax` in [glmnet::glmnet].
#'   this will determine how many non-zero coefficients we would like to have.
#'
#' @return a vector of beta coefficients
#'
#' @noRd
#'
#' @examples
#'
#' penalized_cph(
#'  x_node = as.matrix(pbc_orsf[, c('age', 'bili', 'chol', 'albumin')]),
#'  y_node = as.matrix(pbc_orsf[, c('time', 'status')]),
#'  w_node = rep(1, nrow(pbc_orsf)),
#'  alpha = 1/2,
#'  df_target = 2
#' )

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

 if(is_error(fit)){
  return(matrix(0, nrow=ncol(x_node), ncol=1))
 }

 for(i in seq_along(fit$df)){
  if(fit$df[i] >= df_target || i == length(fit$df)){
   return(matrix(fit$beta[, i, drop=TRUE], ncol = 1))
  }
 }

}

