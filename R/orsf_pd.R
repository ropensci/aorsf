


#' ORSF partial dependence
#'
#' Compute partial dependence for an ORSF model.
#' `r roxy_pd_explain()`
#' `r roxy_pd_oob_explain('partial dependence')`
#'
#' @inheritParams predict.ObliqueForest
#'
#'
#' @param pred_spec (*named list* or _data.frame_).
#'
#'  -  If `pred_spec` is a named list,
#'   Each item in the list should be a vector of values that will be used as
#'   points in the partial dependence function. The name of each item in the
#'   list should indicate which variable will be modified to take the
#'   corresponding values.
#'
#'  - If `pred_spec` is a `data.frame`, columns will
#'   indicate variable names, values will indicate variable values, and
#'   partial dependence will be computed using the inputs on each row.
#'
#' @param pred_type (_character_) the type of predictions to compute. Valid
#'   options are
#'
#'   - 'risk' : probability of having an event at or before `pred_horizon`.
#'   - 'surv' : 1 - risk.
#'   - 'chf': cumulative hazard function
#'   - 'mort': mortality prediction
#'
#' @param na_action `r roxy_na_action_header("new_data")`
#'
#'   - `r roxy_na_action_fail("new_data")`
#'   - `r roxy_na_action_omit("new_data")`
#'
#' @param expand_grid (_logical_) if `TRUE`, partial dependence will be
#'   computed at all possible combinations of inputs in `pred_spec`. If
#'   `FALSE`, partial dependence will be computed for each variable
#'   in `pred_spec`, separately.
#'
#' @param prob_values (_numeric_) a vector of values between 0 and 1,
#'   indicating what quantiles will be used to summarize the partial
#'   dependence values at each set of inputs. `prob_values` should
#'   have the same length as `prob_labels`. The quantiles are calculated
#'   based on predictions from `object` at each set of values indicated
#'   by `pred_spec`.
#'
#' @param prob_labels (_character_) a vector of labels with the same length
#'   as `prob_values`, with each label indicating what the corresponding
#'   value in `prob_values` should be labelled as in summarized outputs.
#'   `prob_labels` should have the same length as `prob_values`.
#'
#' @param boundary_checks (_logical_) if `TRUE`, `pred_spec` will be checked
#'  to make sure the requested values are between the 10th and 90th
#'  percentile in the object's training data. If `FALSE`, these checks are
#'  skipped.
#'
#' @param n_thread `r roxy_n_thread_header("computing predictions")`
#'
#' @param ... `r roxy_dots()`
#'
#' @return a [data.table][data.table::data.table-package] containing
#'   partial dependence values for the specified variable(s) at the
#'   specified prediction horizon(s).
#'
#' @details
#'
#' `r roxy_pd_limitations()`
#'
#' @references
#'
#' `r roxy_cite_hooker_2021()`
#'
#' @includeRmd Rmd/orsf_pd_examples.Rmd
#'
#' @export
#'
orsf_pd_oob <- function(object,
                        pred_spec,
                        pred_horizon = NULL,
                        pred_type = 'risk',
                        expand_grid = TRUE,
                        prob_values = c(0.025, 0.50, 0.975),
                        prob_labels = c('lwr', 'medn', 'upr'),
                        boundary_checks = TRUE,
                        n_thread = 1,
                        ...){

 check_dots(list(...), orsf_pd_oob)

 orsf_pred_dependence(object = object,
                      pred_spec = pred_spec,
                      pd_data = NULL,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      expand_grid = expand_grid,
                      prob_values = prob_values,
                      prob_labels = prob_labels,
                      boundary_checks = boundary_checks,
                      n_thread = n_thread,
                      oobag = TRUE,
                      type_output = 'smry')

}

#' @rdname orsf_pd_oob
#' @export
orsf_pd_inb <- function(object,
                        pred_spec,
                        pred_horizon = NULL,
                        pred_type = 'risk',
                        expand_grid = TRUE,
                        prob_values = c(0.025, 0.50, 0.975),
                        prob_labels = c('lwr', 'medn', 'upr'),
                        boundary_checks = TRUE,
                        n_thread = 1,
                        ...){

 check_dots(list(...), orsf_pd_inb)

 if(is.null(object$data))
  stop("no data were found in object. ",
       "did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)

 orsf_pred_dependence(object = object,
                      pred_spec = pred_spec,
                      pd_data = object$data,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      expand_grid = expand_grid,
                      prob_values = prob_values,
                      prob_labels = prob_labels,
                      boundary_checks = boundary_checks,
                      n_thread = n_thread,
                      oobag = FALSE,
                      type_output = 'smry')

}

#' @rdname orsf_pd_oob
#' @export
orsf_pd_new <- function(object,
                        pred_spec,
                        new_data,
                        pred_horizon = NULL,
                        pred_type = 'risk',
                        na_action = 'fail',
                        expand_grid = TRUE,
                        prob_values = c(0.025, 0.50, 0.975),
                        prob_labels = c('lwr', 'medn', 'upr'),
                        boundary_checks = TRUE,
                        n_thread = 1,
                        ...){

 check_dots(list(...), orsf_pd_new)

 orsf_pred_dependence(object = object,
                      pred_spec = pred_spec,
                      pd_data = new_data,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      na_action = na_action,
                      expand_grid = expand_grid,
                      prob_values = prob_values,
                      prob_labels = prob_labels,
                      boundary_checks = boundary_checks,
                      n_thread = n_thread,
                      oobag = FALSE,
                      type_output = 'smry')

}


#' ORSF Individual Conditional Expectations
#'
#' Compute individual conditional expectations for an ORSF model.
#' `r roxy_ice_explain()`
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
                         pred_spec,
                         pred_horizon = NULL,
                         pred_type = 'risk',
                         expand_grid = TRUE,
                         boundary_checks = TRUE,
                         n_thread = 1,
                         ...){

 check_dots(list(...), orsf_ice_oob)

 orsf_pred_dependence(object = object,
                      pred_spec = pred_spec,
                      pd_data = NULL,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      expand_grid = expand_grid,
                      boundary_checks = boundary_checks,
                      n_thread = n_thread,
                      oobag = TRUE,
                      type_output = 'ice')

}

#' @rdname orsf_ice_oob
#' @export
orsf_ice_inb <- function(object,
                         pred_spec,
                         pred_horizon = NULL,
                         pred_type = 'risk',
                         expand_grid = TRUE,
                         boundary_checks = TRUE,
                         n_thread = 1,
                         ...){

 check_dots(list(...), orsf_ice_oob)

 if(is.null(object$data))
  stop("no data were found in object. ",
       "did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)

 orsf_pred_dependence(object = object,
                      pred_spec = pred_spec,
                      pd_data = object$data,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      expand_grid = expand_grid,
                      boundary_checks = boundary_checks,
                      n_thread = n_thread,
                      oobag = FALSE,
                      type_output = 'ice')

}

#' @rdname orsf_ice_oob
#' @export
orsf_ice_new <- function(object,
                         pred_spec,
                         new_data,
                         pred_horizon = NULL,
                         pred_type = 'risk',
                         na_action = 'fail',
                         expand_grid = TRUE,
                         boundary_checks = TRUE,
                         n_thread = 1,
                         ...){

 check_dots(list(...), orsf_ice_new)

 orsf_pred_dependence(object = object,
                      pred_spec = pred_spec,
                      pd_data = new_data,
                      pred_horizon = pred_horizon,
                      pred_type = pred_type,
                      na_action = na_action,
                      expand_grid = expand_grid,
                      boundary_checks = boundary_checks,
                      n_thread = n_thread,
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
                                 pred_spec,
                                 pred_horizon,
                                 pred_type,
                                 na_action = 'fail',
                                 expand_grid,
                                 prob_values = NULL,
                                 prob_labels = NULL,
                                 boundary_checks,
                                 n_thread,
                                 oobag,
                                 type_output){

 check_arg_is(object, arg_name = 'object', expected_class = 'ObliqueForest')

 pred_horizon <- infer_pred_horizon(object, pred_type, pred_horizon)

 if(is.null(prob_values)) prob_values <- c(0.025, 0.50, 0.975)
 if(is.null(prob_labels)) prob_labels <- c('lwr', 'medn', 'upr')

 if(oobag && is.null(object$data))
  stop("no data were found in object. ",
       "did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)

 if(is.null(pred_horizon) && pred_type %in% c('risk', 'surv', 'chf')){
  stop("pred_horizon must be specified for ",
       pred_type, " predictions.", call. = FALSE)
 }

 check_arg_type(arg_value = boundary_checks,
                arg_name = 'boundary_checks',
                expected_type = 'logical')

 check_arg_length(arg_value = boundary_checks,
                  arg_name = 'boundary_checks',
                  expected_length = 1)

 check_arg_type(arg_value = expand_grid,
                arg_name = 'expand_grid',
                expected_type = 'logical')

 check_arg_length(arg_value = expand_grid,
                  arg_name = 'expand_grid',
                  expected_length = 1)

 check_arg_type(arg_value = prob_values,
                arg_name = 'prob_values',
                expected_type = 'numeric')

 check_arg_gteq(arg_value = prob_values,
                arg_name = 'prob_values',
                bound = 0)

 check_arg_lteq(arg_value = prob_values,
                arg_name = 'prob_values',
                bound = 1)

 check_arg_type(arg_value = prob_labels,
                arg_name = 'prob_labels',
                expected_type = 'character')

 if(length(prob_values) != length(prob_labels)){
  stop("prob_values and prob_labels must have the same length.",
       call. = FALSE)
 }

 check_arg_type(arg_value = oobag,
                arg_name = 'oobag',
                expected_type = 'logical')

 check_arg_length(arg_value = oobag,
                  arg_name = 'oobag',
                  expected_length = 1)

 check_arg_type(arg_value = pred_type,
                arg_name = "pred_type",
                expected_type = 'character')

 check_arg_length(arg_value = pred_type,
                  arg_name = "pred_type",
                  expected_length = 1)

 check_arg_is_valid(arg_value = pred_type,
                    arg_name = "pred_type",
                    valid_options = c("risk",
                                      "surv",
                                      "chf",
                                      "mort",
                                      "prob"))

 object$compute_dependence(pd_data = pd_data,
                           pred_spec = pred_spec,
                           pred_horizon = pred_horizon,
                           pred_type = pred_type,
                           na_action = na_action,
                           expand_grid = expand_grid,
                           prob_values = prob_values,
                           prob_labels = prob_labels,
                           boundary_checks = boundary_checks,
                           n_thread = n_thread,
                           oobag = oobag,
                           type_output = type_output)


}
