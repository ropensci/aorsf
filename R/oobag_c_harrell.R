

#' Harrell's C-statistic
#'
#' This function is for testing and internal use.
#'
#' @param y_mat outcome matrix
#' @param s_vec vector of predicted survival
#'
#' @return the C-statistic
#'
#' @noRd
#'

oobag_c_survival <- function(y_mat, w_vec, s_vec){

 survival::concordancefit(
  y = survival::Surv(y_mat),
  x = s_vec
 )$concordance

}


