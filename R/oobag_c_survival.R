

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

 data <- as.data.frame(cbind(y_mat, s_vec))
 names(data) = c("time", "status", "x")

 survival::concordance(
  survival::Surv(time, status) ~ x,
  data = data,
  weights = w_vec
 )$concordance

}

oobag_c_risk <- function(y_mat, w_vec, s_vec){

 data <- as.data.frame(cbind(y_mat, s_vec))
 names(data) = c("time", "status", "x")

 1 - survival::concordance(
  survival::Surv(time, status) ~ x,
  data = data,
  weights = w_vec
 )$concordance

}


