


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
                            expand_grid = TRUE,
                            prob_values = c(0.025, 0.50, 0.975),
                            prob_labels = c('lwr', 'est', 'upr'),
                            oobag = TRUE,
                            risk = TRUE){

 # TODO: fast way to do this when going through variables sequentially
 # instead of using expand.grid?

 check_predict(object, pd_data, times, risk)

 x_new <- as.matrix(
  one_hot(x_data = pd_data,
          fi = get_fctr_info(object),
          names_x_data = get_names_x(object))
 )

 pd_smry_fun <- switch(as.character(expand_grid),
                       'TRUE' = orsf_pd_summary_grid,
                       'FALSE' = orsf_pd_summary_loop)

 pd_fun <- switch(paste(length(times) == 1, oobag, sep = "_"),
                  'TRUE_FALSE' = new_pd_smry_uni,
                  'TRUE_TRUE' = oob_pd_smry_uni,
                  'FALSE_FALSE' = stop("Not ready yet"),
                  'FALSE_TRUE' = stop("Not ready yet"))

 pd_smry_fun(object, x_new, pd_spec, times, pd_fun,
             prob_values, prob_labels, oobag, risk)

}


orsf_pd_summary_grid <- function(object, x_new, pd_spec, times,
                                 pd_fun, prob_values, prob_labels,
                                 oobag, risk){

 pd_grid <- expand.grid(pd_spec, stringsAsFactors = TRUE)

 pd_grid_new <- one_hot(x_data = pd_grid,
                        fi = get_fctr_info(object),
                        names_x_data = names(pd_grid))

 x_cols <- match(names(pd_grid_new), colnames(x_new))



 pd_vals <- pd_fun(forest      = object$forest,
                   x_new_      = x_new,
                   x_cols_     = x_cols-1,
                   x_vals_     = as.matrix(pd_grid_new),
                   probs_      = prob_values,
                   time_dbl    = times,
                   return_risk = risk)

 rownames(pd_vals) <- c('mean', prob_labels)

 return(cbind(pd_grid, t(pd_vals)))

}


orsf_pd_summary_loop <- function(object, x_new, pd_spec, times,
                                 pd_fun, prob_values, prob_labels,
                                 oobag, risk){

 fi <- get_fctr_info(object)

 output <- vector(mode = 'list', length = length(pd_spec))

 for(i in seq_along(pd_spec)){

  pd_new  <- as.data.frame(pd_spec[i])
  pd_name <- names(pd_spec)[i]
  pd_bind <- data.frame(name = pd_name, value = as.character(pd_spec[[i]]))

  if(pd_name %in% fi$cols) pd_new <- one_hot(pd_new,
                                             fi = fi,
                                             names_x_data = pd_name)

  x_cols <- match(names(pd_new), colnames(x_new))

  x_vals <- x_new[, x_cols]


  pd_vals <- pd_fun(forest      = object$forest,
                    x_new_      = x_new,
                    x_cols_     = x_cols-1,
                    x_vals_     = as.matrix(pd_new),
                    probs_      = prob_values,
                    time_dbl    = times,
                    return_risk = risk)

  # pd_fun modifies x_new by reference, so reset it.
  x_new[, x_cols] <- x_vals

  rownames(pd_vals) <- c('mean', prob_labels)

  output[[i]] <- cbind(pd_bind, t(pd_vals))


 }

 Reduce(rbind, output)

}


