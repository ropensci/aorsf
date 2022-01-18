


#' ORSF partial dependence
#'
#' @srrstats {G1.4} *documented with Roxygen*
#'
#' @inheritParams predict.aorsf
#'
#' @param object (_aorsf_) An accelerated oblique random survival forest model.
#'
#' @param pd_data (_data frame_) that will be used to compute partial
#'   dependence. If `NULL`, then the training data of `object` will be
#'   used. If the training data were not attached to `object`
#'   (see `attach_data` input in [orsf]), an error will be triggered.
#'
#' @param pd_spec (_named list_ or _data.frame_). If `pd_spec` is a named list,
#'   Each item in the list should be a vector of values that will be used as
#'   points in the partial dependence function. The name of each item in the
#'   list should indicate which variable will be modified to take the
#'   corresponding values. If `pd_spec` is a `data.frame`, columns will
#'   indicate variable names, values will indicate variable values, and
#'   partial dependence will be computed using the inputs on each row.
#'
#' @param expand_grid (_logical_) if `TRUE`, partial dependence will be
#'   computed at all possible combinations of inputs in `pd_spec`. If
#'   `FALSE`, partial dependence will be computed for each variable
#'   in `pd_spec`, separately.
#'
#' @param prob_values (_numeric_) a vector of values between 0 and 1,
#'   indicating what quantiles will be used to summarize the partial
#'   dependence values at each set of inputs.
#'
#' @param prob_labels (_character_) a vector of labels with the same length
#'   as `prob_values`, with each label indicating what the corresponding
#'   value in `prob_values` should be labelled as in summarized outputs.
#'
#' @param oobag (_logical_) if `TRUE`, then partial dependence will be
#'   computed using the out of bag training data. You should set
#'   `oobag = TRUE` if you are computing partial dependence using the
#'   training data for `object`.
#'
#' @param boundary_checks (_logical_) if `TRUE`, `pd_spec` will be vetted
#'  to make sure the requested values are between the 10th and 90th
#'  percentile in the object's training data. If `FALSE`, these checks are
#'  skipped.
#'
#' @return a `data.table` containing summarized partial dependence
#'   values if using `orsf_pd_summery` or individual conditional
#'   expectation (ICE) partial dependence if using `orsf_pd_ice`.
#'
#' @export
#'
#' @examples
#'
#' fit <- orsf(pbc_orsf, Surv(time, status) ~ . - id)
#'
#' orsf_pd_summary(fit, pd_spec = list(bili = c(1,2,3,4,5,6)), times = 1000)
#'
#' # more points for a plot
#' pd_spec <- list(bili = seq(1, 6, length.out = 20))
#' data_ice <- orsf_pd_ice(fit, pd_spec = pd_spec, times = 1000)
#'
#' head(data_ice)
#'
#' library(ggplot2)
#' ggplot(data_ice) +
#'  aes(x = bili, y = pred, group = id_row) +
#'  geom_line(alpha = 0.4, color = 'grey') +
#'  geom_smooth(aes(group = 1), color = 'black', se = FALSE) +
#'  theme(panel.grid = element_blank())

orsf_pd_summary <- function(object,
                            pd_data = NULL,
                            pd_spec,
                            times = NULL,
                            expand_grid = TRUE,
                            prob_values = c(0.025, 0.50, 0.975),
                            prob_labels = c('lwr', 'medn', 'upr'),
                            oobag = TRUE,
                            risk = TRUE,
                            boundary_checks = TRUE){


 check_pd_inputs(object = object,
                 expand_grid = expand_grid,
                 prob_values = prob_values,
                 prob_labels = prob_labels,
                 oobag = oobag,
                 risk = risk)

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
          risk = risk,
          boundary_checks = boundary_checks)

}

#' @rdname orsf_pd_summary
#' @export
orsf_pd_ice <- function(object,
                        pd_data = NULL,
                        pd_spec,
                        times = NULL,
                        expand_grid = TRUE,
                        oobag = TRUE,
                        risk = TRUE,
                        boundary_checks = TRUE){


 check_pd_inputs(object = object,
                 expand_grid = expand_grid,
                 oobag = oobag,
                 risk = risk)

 orsf_pd_(object = object,
          pd_data = pd_data,
          pd_spec = pd_spec,
          times = times,
          type_output = 'ice',
          type_input = if(expand_grid) 'grid' else 'loop',
          prob_values = c(0.025, 0.50, 0.975),
          prob_labels = c('lwr', 'medn', 'upr'),
          oobag = oobag,
          risk = risk,
          boundary_checks = boundary_checks)

}


#' delegation function of orsf_pd family
#'
#' this function takes inputs from the main API functions and
#'   determines which of the lower level functions to call.
#'
#' @inheritParams orsf_pd_summary
#' @param type_output 'ice' or 'smry'.
#'   this in combination with oobag determines which cpp routine to use.
#' @param type_input if 'grid', then all combos of pd_data are considered.
#'   if 'loop', then each entry of pd_data is considered separately.
#'
#' @return output from one of the pd working functions
#'
#' @noRd

orsf_pd_ <- function(object,
                     pd_data,
                     pd_spec,
                     times,
                     type_output,
                     type_input,
                     prob_values,
                     prob_labels,
                     oobag,
                     risk,
                     boundary_checks){

 if(is.null(times)) times <- object$time_pred

 if(is.null(times))
  stop("times was not specified and could not be found in object. ",
       "did you use oobag_pred = FALSE when running orsf()?",
       call. = FALSE)

 if(length(times) > 1){
  stop("orsf_pd functions only allow 1 prediction time,",
       " but your times input has length ", length(times), ".",
       call. = FALSE)
 }

 if(oobag && !is.null(pd_data)){

  warning(
   "pd_data should be NULL when computing out-of-bag partial dependence,",
   " i.e., when oobag = TRUE. the pd_data input will be ignored and instead",
   " object$data_train will be used to compute out-of-bag partial dependence.",
   " If you want to compute partial dependence on new data, set oobag = FALSE",
   call. = FALSE
  )

  pd_data <- object$data_train

 }



 if(is.null(pd_data)) pd_data <- object$data_train

 if(is.null(pd_data)) stop("training data were not found in object. ",
                           "did you use attach_data = FALSE when ",
                           "running orsf()?", call. = FALSE)



 check_predict(object, pd_data, times, risk)

 if(is_empty(pd_spec)){

   stop("pd_spec is empty", call. = FALSE)

 }

 if(is_empty(names(pd_spec))){

  stop("pd_spec is unnamed", call. = FALSE)

 }

 bad_name_index <- which(is.na(match(names(pd_spec), get_names_x(object))))

 if(!is_empty(bad_name_index)){

  bad_names <- names(pd_spec)[bad_name_index]

  stop("some variables in pd_spec are not in object's training data: ",
       paste_collapse(bad_names, last = ' and '),
       call. = FALSE)

 }

 numeric_bounds <- get_numeric_bounds(object)
 numeric_names <- intersect(colnames(numeric_bounds), names(pd_spec))

 if(!is_empty(numeric_names) && boundary_checks){

  for(.name in numeric_names){

   vals_above_stop <- which(pd_spec[[.name]] > numeric_bounds['90%', .name])
   vals_below_stop <- which(pd_spec[[.name]] < numeric_bounds['10%', .name])

   boundary_error <- FALSE
   vals_above_list <- vals_below_list <- " "

   if(!is_empty(vals_above_stop)){
    vals_above_list <- paste_collapse(
     table.glue::table_value(pd_spec[[.name]][vals_above_stop]),
     last = ' and '
    )

    boundary_error <- TRUE
    vals_above_list <-
     paste0(" (",vals_above_list," > ", numeric_bounds['90%', .name],") ")

   }

   if(!is_empty(vals_below_stop)){

    vals_below_list <- paste_collapse(
     table.glue::table_value(pd_spec[[.name]][vals_below_stop]),
     last = ' and '
    )

    boundary_error <- TRUE

    vals_below_list <-
       paste0(" (",vals_below_list," < ", numeric_bounds['10%', .name],") ")

   }

   if(boundary_error)
    stop("Some values for ",
         .name,
         " in pd_spec are above",
         vals_above_list,
         "or below",
         vals_below_list,
         "90th or 10th percentiles in training data.",
         " Change pd_spec or set boundary_checks = FALSE",
         " to prevent this error",
         call. = FALSE)

  }

 }


 x_new <- as.matrix(
  ref_code(x_data = pd_data,
          fi = get_fctr_info(object),
          names_x_data = get_names_x(object))
 )


 if(is.data.frame(pd_spec)) type_input <- 'grid'

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


#' grid working function in orsf_pd family
#'
#' This function expands pd_spec into a grid with all combos of inputs,
#'   and computes partial dependence for each one.
#'
#' @inheritParams orsf_pd_
#' @param x_new the x-matrix used to compute partial dependence
#' @param pd_fun_predict which cpp function to use.
#'
#' @return a `data.table` containing summarized partial dependence
#'   values if using `orsf_pd_summery` or individual conditional
#'   expectation (ICE) partial dependence if using `orsf_pd_ice`.
#'
#' @noRd

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

 if(!is.data.frame(pd_spec))
  pd_spec <- expand.grid(pd_spec, stringsAsFactors = TRUE)

 fi_ref <- get_fctr_info(object)

 for(i in seq_along(fi_ref$cols)){

  ii <- fi_ref$cols[i]

  if(is.character(pd_spec[[ii]]) && !fi_ref$ordr[i]){

   pd_spec[[ii]] <- factor(pd_spec[[ii]],
                           levels = fi_ref$lvls[[ii]])

  }

 }

 check_new_data_fctrs(new_data  = pd_spec,
                      names_x   = get_names_x(object),
                      fi_ref    = fi_ref,
                      label_new = "pd_spec")

 pd_spec_new <- ref_code(x_data = pd_spec,
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
  .names <- names(output)

 }

 if(type_output == 'ice'){

  colnames(pd_vals) <- c('id_variable', 'pred')
  pd_spec$id_variable <- seq(nrow(pd_spec))
  output <- merge(pd_spec, pd_vals, by = 'id_variable')
  output$id_row <- rep(seq(nrow(x_new)), times = nrow(pd_spec))

  ids <- c('id_variable', 'id_row')
  .names <- c(ids, setdiff(names(output), ids))

 }

 as.data.table(output[, .names])

}

#' loop working function in orsf_pd family
#'
#' This function loops through the items in pd_spec one by one,
#'   computing partial dependence for each one separately.
#'
#' @inheritParams orsf_pd_
#' @param x_new the x-matrix used to compute partial dependence
#' @param pd_fun_predict which cpp function to use.
#'
#' @return a `data.table` containing summarized partial dependence
#'   values if using `orsf_pd_summery` or individual conditional
#'   expectation (ICE) partial dependence if using `orsf_pd_ice`.
#'
#' @noRd

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

  pd_bind <- data.frame(variable = pd_name,
                        value = rep(NA_real_, length(pd_spec[[i]])),
                        level = rep(NA_character_, length(pd_spec[[i]])))

  if(pd_name %in% fi$cols) {

   pd_bind$level <- as.character(pd_spec[[i]])

   pd_new <- ref_code(pd_new,
                     fi = fi,
                     names_x_data = pd_name)
  } else {

   pd_bind$value <- pd_spec[[i]]

  }

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

   colnames(pd_vals) <- c('id_variable', 'pred')
   pd_bind$id_variable <- seq(nrow(pd_bind))
   output[[i]] <- merge(pd_bind, pd_vals, by = 'id_variable')
   output[[i]]$id_row <- seq(nrow(output[[i]]))

  }

 }

 output <- rbindlist(output)

 if(type_output == 'ice'){

  ids <- c('id_variable', 'id_row')
  .names <- c(ids, setdiff(names(output), ids))
  setcolorder(output, neworder = .names)

 }


 output

}

