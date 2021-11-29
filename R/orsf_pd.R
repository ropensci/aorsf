


#' Title
#'
#' @param object
#' @param pd_data
#' @param pd_spec
#' @param times
#' @param probs
#' @param risk
#'
#' @return
#' @export
#'
#' @examples
orsf_pd_summary <- function(object,
                            pd_data,
                            pd_spec,
                            times,
                            probs = c(0.025, 0.975),
                            risk = TRUE){

 check_predict(object, pd_data, times, risk)

 pd_grid <- expand.grid(pd_spec)

 x_new <- as.matrix(
  one_hot(data = new_data,
          fi = get_fctr_info(object),
          names_x_data = get_names_x(object))
 )

 # pd_mean <- pd_lwr <- pd_upr <- rep(NA_real_, nrow(pd_grid))
 #
 # for(i in seq(nrow(pd_grid))){
 #
 #  for(j in seq(ncol(pd_grid))){
 #   pd_data[, names(pd_grid)[j]] <- pd_grid[i, j]
 #  }
 #
 #  pd_vals <- predict(object,
 #                     new_data = pd_data,
 #                     times = times,
 #                     risk = risk)
 #
 #  pd_quant <- quantile(pd_vals, probs = probs)
 #
 #  pd_mean[i] <- mean(pd_vals)
 #  pd_lwr[i] <- pd_quant[1]
 #  pd_upr[i] <- pd_quant[2]
 #
 # }
 #
 # cbind(pd_grid, mean = pd_mean, lwr = pd_lwr, upr = pd_upr)

}
