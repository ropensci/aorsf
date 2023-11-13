

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
#'
#' penalized_logreg(
#'  x_node = as.matrix(penguins_orsf[, c('bill_length_mm',
#'                                       'bill_depth_mm',
#'                                       'flipper_length_mm',
#'                                       'body_mass_g')]),
#'  y_node = as.matrix(as.numeric(penguins_orsf$species == 'Adelie')),
#'  w_node = rep(1, nrow(penguins_orsf)),
#'  alpha = 1/2,
#'  df_target = 2
#' )



penalized_cph <- function(x_node,
                          y_node,
                          w_node,
                          alpha,
                          df_target){

 colnames(y_node) <- c('time', 'status')

 penalized_fitter(x_node = x_node,
                  y_node = y_node,
                  w_node = w_node,
                  alpha = alpha,
                  df_target = df_target,
                  family = "cox")

}

penalized_logreg <- function(x_node,
                             y_node,
                             w_node,
                             alpha,
                             df_target){

 col <- sample(ncol(y_node), 1)
 y_node <- as.factor(y_node[, col])
 w_node <- as.numeric(w_node)

 penalized_fitter(x_node = x_node,
                  y_node = y_node,
                  w_node = w_node,
                  alpha = alpha,
                  df_target = df_target,
                  family = "binomial")

}

penalized_linreg <- function(x_node,
                             y_node,
                             w_node,
                             alpha,
                             df_target){

 w_node <- as.numeric(w_node)

 penalized_fitter(x_node = x_node,
                  y_node = y_node,
                  w_node = w_node,
                  alpha = alpha,
                  df_target = df_target,
                  family = "gaussian")

}

penalized_fitter <- function(x_node,
                             y_node,
                             w_node,
                             alpha,
                             df_target,
                             family){

 suppressWarnings(
  fit <- try(
   glmnet::glmnet(x = x_node,
                  y = y_node,
                  weights = w_node,
                  alpha = alpha,
                  family = family),
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

