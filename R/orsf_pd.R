


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
                            prob_values = c(0.025, 0.50, 0.975),
                            prob_labels = c('lwr', 'est', 'upr'),
                            risk = TRUE){

 check_predict(object, pd_data, times, risk)

 pd_grid <- expand.grid(pd_spec)

 x_new <- as.matrix(
  one_hot(x_data = pd_data,
          fi = get_fctr_info(object),
          names_x_data = get_names_x(object))
 )

 pd_grid_new <- one_hot(x_data = pd_grid,
                        fi = get_fctr_info(object),
                        names_x_data = names(pd_grid))

 x_cols <- match(names(pd_grid_new), colnames(x_new))

 if(length(times) == 1){

  pd_vals <- orsf_pd_smry_uni(forest      = object$forest,
                              x_new_      = x_new,
                              x_cols_     = x_cols-1,
                              x_vals_     = as.matrix(pd_grid_new),
                              probs_      = prob_values,
                              time_dbl    = times,
                              return_risk = risk)

  rownames(pd_vals) <- c('mean', prob_labels)

  return(cbind(pd_grid, t(pd_vals)))

 }

}
