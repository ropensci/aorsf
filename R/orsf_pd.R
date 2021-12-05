


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
                            pd_data = NULL,
                            pd_spec,
                            times,
                            expand_grid = TRUE,
                            prob_values = c(0.025, 0.50, 0.975),
                            prob_labels = c('lwr', 'est', 'upr'),
                            oobag = TRUE,
                            risk = TRUE){

 check_call(
  match.call(),
  expected = list(
   object = list(
    class = 'aorsf'
   ),
   'expand_grid' = list(
    type = 'logical',
    length = 1
   ),
   prob_values = list(
    type = 'numeric',
    lwr = 0,
    upr = 1
   ),
   prob_labels = list(
    type = 'character'
   )
  )
 )

 if(length(prob_values) != length(prob_labels)){
  stop("prob_values and prob_labels must have the same length.",
       call. = FALSE)
 }

 orsf_pd_(object = object,
          pd_data = pd_data,
          pd_spec = pd_spec,
          times = times,
          type_output = 'smry',
          type_input = if(expand_grid) 'grid' else 'loop',
          prob_values = prob_values,
          prob_labels = prob_labels,
          oobag = oobag,
          risk = risk)

}

#' @rdname orsf_pd_summary
#' @export
orsf_pd_ice <- function(object,
                        pd_data = NULL,
                        pd_spec,
                        times,
                        expand_grid = TRUE,
                        oobag = TRUE,
                        risk = TRUE){

 check_call(
  match.call(),
  expected = list(
   object = list(
    class = 'aorsf'
   ),
   'expand_grid' = list(
    type = 'logical',
    length = 1
   )
  )
 )

 orsf_pd_(object = object,
          pd_data = pd_data,
          pd_spec = pd_spec,
          times = times,
          type_output = 'ice',
          type_input = if(expand_grid) 'grid' else 'loop',
          prob_values = prob_values,
          prob_labels = prob_labels,
          oobag = oobag,
          risk = risk)

}


orsf_pd_ <- function(object,
                     pd_data,
                     pd_spec,
                     times,
                     type_output,
                     type_input,
                     prob_values,
                     prob_labels,
                     oobag,
                     risk){


 if(length(times) > 1){
  stop("orsf_pd functions only allow 1 prediction time,",
       " but your times input has length ", length(times), ".",
       call. = FALSE)
 }

 if(is.null(pd_data)) pd_data <- object$data_train

 if(is.null(pd_data)) stop("training data were not found in object. ",
                           "did you use attach_data = FALSE when ",
                           "running orsf()?", call. = FALSE)

 check_predict(object, pd_data, times, risk)

 Call <- match.call()

 check_call(
  Call,
  expected = list(
   'prob_values' = list(
    type = 'numeric',
    lwr = 0,
    upr = 1
   ),
   'prob_labels' = list(
    type = 'character'
   ),
   'oobag' = list(
    type = 'logical',
    length = 1
   ),
   'risk' = list(
    type = 'logical',
    length = 1
   )
  )
 )

 if(is_empty(pd_spec)){

   stop("pd_spec is empty", call. = FALSE)

 }

 bad_name_index <- which(is.na(match(names(pd_spec), get_names_x(object))))

 if(!is_empty(bad_name_index)){

  bad_names <- names(pd_spec)[bad_name_index]

  stop("some variables in pd_spec are not in object's training data: ",
       paste_collapse(bad_names, last = ' and '),
       call. = FALSE)

 }


 x_new <- as.matrix(
  one_hot(x_data = pd_data,
          fi = get_fctr_info(object),
          names_x_data = get_names_x(object))
 )


 pd_fun_structure <- switch(type_input,
                            'grid' = pd_grid,
                            'loop' = pd_loop)

 pd_fun_predict <- switch(paste(oobag, type_output, sep = "_"),
                          "TRUE_ice" = pd_oob_ice,
                          "TRUE_smry" = pd_oob_smry,
                          "FALSE_ice" = pd_new_ice,
                          "FALSE_smry" = pd_new_smry)


 pd_fun_structure(object,
                  x_new,
                  pd_spec,
                  times,
                  pd_fun_predict,
                  type_output,
                  prob_values,
                  prob_labels,
                  oobag,
                  risk)

}


pd_grid <- function(object,
                    x_new,
                    pd_spec,
                    times,
                    pd_fun_predict,
                    type_output,
                    prob_values,
                    prob_labels,
                    oobag,
                    risk){

 pd_spec <- expand.grid(pd_spec, stringsAsFactors = TRUE)

 pd_spec_new <- one_hot(x_data = pd_spec,
                        fi = get_fctr_info(object),
                        names_x_data = names(pd_spec))

 x_cols <- match(names(pd_spec_new), colnames(x_new))

 pd_vals <- pd_fun_predict(forest      = object$forest,
                           x_new_      = x_new,
                           x_cols_     = x_cols-1,
                           x_vals_     = as.matrix(pd_spec_new),
                           probs_      = prob_values,
                           time_dbl    = times,
                           return_risk = risk)


 if(type_output == 'smry'){

  rownames(pd_vals) <- c('mean', prob_labels)
  output <- cbind(pd_spec, t(pd_vals))

 }

 if(type_output == 'ice'){

  colnames(pd_vals) <- c('key', 'pred')
  pd_spec$key <- seq(nrow(pd_spec))
  output <- merge(pd_spec, pd_vals, by = 'key')
  output$key <- NULL

 }

 output

}


pd_loop <- function(object,
                    x_new,
                    pd_spec,
                    times,
                    pd_fun_predict,
                    type_output,
                    prob_values,
                    prob_labels,
                    oobag,
                    risk){

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


  pd_vals <- pd_fun_predict(forest      = object$forest,
                            x_new_      = x_new,
                            x_cols_     = x_cols-1,
                            x_vals_     = as.matrix(pd_new),
                            probs_      = prob_values,
                            time_dbl    = times,
                            return_risk = risk)

  # pd_fun_predict modifies x_new by reference, so reset it.
  x_new[, x_cols] <- x_vals

  if(type_output == 'smry'){

   rownames(pd_vals) <- c('mean', prob_labels)
   output[[i]] <- cbind(pd_bind, t(pd_vals))

  }

  if(type_output == 'ice'){

   colnames(pd_vals) <- c('key', 'pred')
   pd_bind$key <- seq(nrow(pd_bind))
   output[[i]] <- merge(pd_bind, pd_vals, by = 'key')
   output[[i]]$key <- NULL

  }

 }

 Reduce(rbind, output)

}

# pd_ice_grid <- function(object, x_new, pd_spec, times,
#                         pd_fun_predict, prob_values, prob_labels,
#                         oobag, risk){
#
#  pd_grid <- expand.grid(pd_spec, stringsAsFactors = TRUE)
#
#  pd_grid_new <- one_hot(x_data = pd_grid,
#                         fi = get_fctr_info(object),
#                         names_x_data = names(pd_grid))
#
#  x_cols <- match(names(pd_grid_new), colnames(x_new))
#
#  pd_vals <- pd_fun_predict(forest      = object$forest,
#                            x_new_      = x_new,
#                            x_cols_     = x_cols-1,
#                            x_vals_     = as.matrix(pd_grid_new),
#                            probs_      = prob_values,
#                            time_dbl    = times,
#                            return_risk = risk)
#
#  browser()
#
#  rownames(pd_vals) <- c('mean', prob_labels)
#
#  return(cbind(pd_grid, t(pd_vals)))
#
# }
#
#
# pd_ice_loop <- function(object, x_new, pd_spec, times,
#                         pd_fun_predict, prob_values, prob_labels,
#                         oobag, risk){
#
#  fi <- get_fctr_info(object)
#
#  output <- vector(mode = 'list', length = length(pd_spec))
#
#  for(i in seq_along(pd_spec)){
#
#   pd_new  <- as.data.frame(pd_spec[i])
#   pd_name <- names(pd_spec)[i]
#   pd_bind <- data.frame(name = pd_name, value = as.character(pd_spec[[i]]))
#
#   if(pd_name %in% fi$cols) pd_new <- one_hot(pd_new,
#                                              fi = fi,
#                                              names_x_data = pd_name)
#
#   x_cols <- match(names(pd_new), colnames(x_new))
#
#   x_vals <- x_new[, x_cols]
#
#
#   pd_vals <- pd_fun_predict(forest      = object$forest,
#                             x_new_      = x_new,
#                             x_cols_     = x_cols-1,
#                             x_vals_     = as.matrix(pd_new),
#                             probs_      = prob_values,
#                             time_dbl    = times,
#                             return_risk = risk)
#
#   browser()
#
#   # pd_fun_predict modifies x_new by reference, so reset it.
#   x_new[, x_cols] <- x_vals
#
#   rownames(pd_vals) <- c('mean', prob_labels)
#
#   output[[i]] <- cbind(pd_bind, t(pd_vals))
#
#
#  }
#
#  Reduce(rbind, output)
#
# }


