


#' ORSF partial dependence
#'
#' @srrstats {G2.1a} *explicit secondary documentation of expectations on data types of all vector inputs*
#' @srrstats {G2.0a} *documenting length by indicating these inputs are vectors and that they must have the same length as the other.*
#'
#' Compute partial dependence for an ORSF model.
#' `r roxy_pd_explain()`
#'
#' `r roxy_pd_oob_explain('partial dependence')`
#'
#' @inheritParams predict.orsf_fit
#'
#' @param pd_spec (*named list* or _data.frame_).
#'
#'  -  If `pd_spec` is a named list,
#'   Each item in the list should be a vector of values that will be used as
#'   points in the partial dependence function. The name of each item in the
#'   list should indicate which variable will be modified to take the
#'   corresponding values.
#'
#'  - If `pd_spec` is a `data.frame`, columns will
#'   indicate variable names, values will indicate variable values, and
#'   partial dependence will be computed using the inputs on each row.
#'
#' @param expand_grid (_logical_) if `TRUE`, partial dependence will be
#'   computed at all possible combinations of inputs in `pd_spec`. If
#'   `FALSE`, partial dependence will be computed for each variable
#'   in `pd_spec`, separately.
#'
#'
#' @param prob_values (_numeric_) a vector of values between 0 and 1,
#'   indicating what quantiles will be used to summarize the partial
#'   dependence values at each set of inputs. `prob_values` should
#'   have the same length as `prob_labels`. The quantiles are calculated
#'   based on predictions from `object` at each set of values indicated
#'   by `pd_spec`.
#'
#' @param prob_labels (_character_) a vector of labels with the same length
#'   as `prob_values`, with each label indicating what the corresponding
#'   value in `prob_values` should be labelled as in summarized outputs.
#'   `prob_labels` should have the same length as `prob_values`.
#'
#' @param boundary_checks (_logical_) if `TRUE`, `pd_spec` will be vetted
#'  to make sure the requested values are between the 10th and 90th
#'  percentile in the object's training data. If `FALSE`, these checks are
#'  skipped.
#'
#' @param ... `r roxy_dots()`
#'
#' @return a [data.table][data.table::data.table-package] containing
#'   partial dependence values for the specified variable(s) at the
#'   specified prediction horizon(s).
#'
#' @includeRmd Rmd/orsf_pd_examples.Rmd
#'
#' @export
#'
#'
orsf_pd_oob <- function(object,
                        pd_spec,
                        pred_horizon = NULL,
                        pred_type = 'risk',
                        expand_grid = TRUE,
                        prob_values = c(0.025, 0.50, 0.975),
                        prob_labels = c('lwr', 'medn', 'upr'),
                        boundary_checks = TRUE,
                        ...){

 check_dots(list(...), orsf_pd_oob)

 orsf_pred_dependence(object = object,
                      pd_spec = pd_spec,
                      pd_data = NULL,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      expand_grid = expand_grid,
                      prob_values = prob_values,
                      prob_labels = prob_labels,
                      boundary_checks = boundary_checks,
                      oobag = TRUE,
                      type_output = 'smry')

}

#' @rdname orsf_pd_oob
#' @export
orsf_pd_inb <- function(object,
                        pd_spec,
                        pred_horizon = NULL,
                        pred_type = 'risk',
                        expand_grid = TRUE,
                        prob_values = c(0.025, 0.50, 0.975),
                        prob_labels = c('lwr', 'medn', 'upr'),
                        boundary_checks = TRUE,
                        ...){

 check_dots(list(...), orsf_pd_inb)

 if(is.null(object$data))
  stop("no data were found in object. ",
       "did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)

 orsf_pred_dependence(object = object,
                      pd_spec = pd_spec,
                      pd_data = object$data,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      expand_grid = expand_grid,
                      prob_values = prob_values,
                      prob_labels = prob_labels,
                      boundary_checks = boundary_checks,
                      oobag = FALSE,
                      type_output = 'smry')

}

#' @rdname orsf_pd_oob
#' @export
orsf_pd_new <- function(object,
                        pd_spec,
                        new_data,
                        pred_horizon = NULL,
                        pred_type = 'risk',
                        expand_grid = TRUE,
                        prob_values = c(0.025, 0.50, 0.975),
                        prob_labels = c('lwr', 'medn', 'upr'),
                        boundary_checks = TRUE,
                        ...){

 check_dots(list(...), orsf_pd_new)

 orsf_pred_dependence(object = object,
                      pd_spec = pd_spec,
                      pd_data = new_data,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      expand_grid = expand_grid,
                      prob_values = prob_values,
                      prob_labels = prob_labels,
                      boundary_checks = boundary_checks,
                      oobag = FALSE,
                      type_output = 'smry')

}


#' ORSF Individual Conditional Expectations
#'
#' Compute individual conditional expectations for an ORSF model.
#' `r roxy_ice_explain()`
#'
#' `r roxy_pd_oob_explain('individual conditional expectations')`
#'
#' @inheritParams orsf_pd_oob
#'
#' @return a [data.table][data.table::data.table-package] containing
#'   individual conditional expectations for the specified variable(s) at the
#'   specified prediction horizon(s).
#'
#' @includeRmd Rmd/orsf_ice_examples.Rmd
#'
#' @export
#'
#'
orsf_ice_oob <- function(object,
                         pd_spec,
                         pred_horizon = NULL,
                         pred_type = 'risk',
                         expand_grid = TRUE,
                         boundary_checks = TRUE,
                         ...){

 check_dots(list(...), orsf_ice_oob)

 orsf_pred_dependence(object = object,
                      pd_spec = pd_spec,
                      pd_data = NULL,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      expand_grid = expand_grid,
                      boundary_checks = boundary_checks,
                      oobag = TRUE,
                      type_output = 'ice')

}

#' @rdname orsf_ice_oob
#' @export
orsf_ice_inb <- function(object,
                         pd_spec,
                         pred_horizon = NULL,
                         pred_type = 'risk',
                         expand_grid = TRUE,
                         boundary_checks = TRUE,
                         ...){

 check_dots(list(...), orsf_ice_oob)

 if(is.null(object$data))
  stop("no data were found in object. ",
       "did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)

 orsf_pred_dependence(object = object,
                      pd_spec = pd_spec,
                      pd_data = object$data,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      expand_grid = expand_grid,
                      boundary_checks = boundary_checks,
                      oobag = FALSE,
                      type_output = 'ice')

}

#' @rdname orsf_ice_oob
#' @export
orsf_ice_new <- function(object,
                         pd_spec,
                         new_data,
                         pred_horizon = NULL,
                         pred_type = 'risk',
                         expand_grid = TRUE,
                         boundary_checks = TRUE,
                         ...){

 check_dots(list(...), orsf_ice_new)

 orsf_pred_dependence(object = object,
                      pd_spec = pd_spec,
                      pd_data = new_data,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      expand_grid = expand_grid,
                      boundary_checks = boundary_checks,
                      oobag = FALSE,
                      type_output = 'ice')


}


#' delegates the work functions for partial dependence
#'
#' this function takes inputs from the main API functions and
#'   determines which of the lower level functions to call.
#'
#' @inheritParams orsf_pd
#' @param type_output 'ice' or 'smry'.
#'   this in combination with oobag determines which cpp routine to use.
#' @param type_input if 'grid', then all combos of pd_data are considered.
#'   if 'loop', then each entry of pd_data is considered separately.
#'
#' @return output from one of the pd working functions
#'
#' @noRd

orsf_pred_dependence <- function(object,
                                 pd_data,
                                 pd_spec,
                                 pred_horizon,
                                 pred_type,
                                 expand_grid,
                                 prob_values = NULL,
                                 prob_labels = NULL,
                                 oobag,
                                 type_output,
                                 boundary_checks){

 pred_horizon <- infer_pred_horizon(object, pred_horizon)

 if(is.null(prob_values)) prob_values <- c(0.025, 0.50, 0.975)
 if(is.null(prob_labels)) prob_labels <- c('lwr', 'medn', 'upr')

 check_pd_inputs(object          = object,
                 pd_spec         = pd_spec,
                 expand_grid     = expand_grid,
                 prob_values     = prob_values,
                 prob_labels     = prob_labels,
                 boundary_checks = boundary_checks,
                 new_data        = pd_data,
                 pred_type       = pred_type,
                 pred_horizon    = pred_horizon)

 if(oobag && is.null(object$data))
  stop("no data were found in object. ",
       "did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)

 if(oobag) pd_data <- object$data

 type_input <- if(expand_grid) 'grid' else 'loop'

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

 risk <- pred_type == 'risk'

 out_list <- lapply(

  X = pred_horizon,

  FUN = function(.pred_horizon){

   pd_fun_structure(object,
                    x_new,
                    pd_spec,
                    .pred_horizon,
                    pd_fun_predict,
                    type_output,
                    prob_values,
                    prob_labels,
                    oobag,
                    risk)

  }

 )

 names(out_list) <- as.character(pred_horizon)

 out <- rbindlist(l = out_list,
                  fill = TRUE,
                  idcol = 'pred_horizon')

 out[, pred_horizon := as.numeric(pred_horizon)]

 # silent print after modify in place
 out[]

 out

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
#'   expectation (ICE) partial dependence if using `orsf_ice`.
#'
#' @noRd

pd_grid <- function(object,
                    x_new,
                    pd_spec,
                    pred_horizon,
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
                           time_dbl    = pred_horizon,
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
  output$id_row <- rep(seq(nrow(x_new)), pred_horizon = nrow(pd_spec))

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
#'   expectation (ICE) partial dependence if using `orsf_ice`.
#'
#' @noRd

pd_loop <- function(object,
                    x_new,
                    pd_spec,
                    pred_horizon,
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
                            time_dbl    = pred_horizon,
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

