

#' Automatic variable values for dependence
#'
#' For partial dependence and individual conditional expectations,
#'   this function allows a variable to be considered without having
#'   to specify what values to set the variable at. The values used
#'   are based on quantiles for continuous variables (10th, 25th, 50th,
#'   75th, and 90th) and unique categories for categorical variables.
#'
#' @param ... names of the variables to use. These can be in quotes
#'   or not in quotes (see examples).
#'
#' @return a character vector with the names
#'
#' @details This function should only be used in the context of
#'   `orsf_pd` or `orsf_ice` functions.
#'
#'
#' @export
#'
#' @examples
#'
#' fit <- orsf(penguins_orsf, species ~., n_tree = 5)
#'
#' orsf_pd_oob(fit, pred_spec_auto(flipper_length_mm))
#'
pred_spec_auto <- function(...){

 input_string <- gsub(".*\\((.*)\\).*", "\\1", match.call())[-1]
 result <- trimws(unlist(strsplit(input_string, ",")))
 class(result) <- c("pspec_auto", class(result))

 result

}


#' Partial dependence
#'
#' Compute partial dependence for an oblique random forest.
#' `r roxy_pd_explain()`
#' `r roxy_pd_oob_explain('partial dependence')`
#'
#' @inheritParams predict.ObliqueForest
#'
#'
#' @param pred_spec (*named list*, *pspec_auto*, or *data.frame*).
#'
#'  -  If `pred_spec` is a named list,
#'   Each item in the list should be a vector of values that will be used as
#'   points in the partial dependence function. The name of each item in the
#'   list should indicate which variable will be modified to take the
#'   corresponding values.
#'
#'  - If `pred_spec` is created using `pred_spec_auto()`, all that is needed
#'    is the names of variables to use (see [pred_spec_auto]).
#'
#'  - If `pred_spec` is a `data.frame`, columns will
#'   indicate variable names, values will indicate variable values, and
#'   partial dependence will be computed using the inputs on each row.
#'
#' @param pred_type (_character_) the type of predictions to compute. Valid
#'   Valid options for survival are:
#'
#'   - 'risk' : probability of having an event at or before `pred_horizon`.
#'   - 'surv' : 1 - risk.
#'   - 'chf': cumulative hazard function
#'   - 'mort': mortality prediction
#'
#'  For classification:
#'
#'  - 'prob': probability for each class
#'
#'  For regression:
#'
#'  - 'mean': predicted mean, i.e., the expected value
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
#' @param verbose_progress (_logical_) if `TRUE`, progress will be
#'  printed to console. If `FALSE` (the default), nothing will be
#'  printed.
#'
#' @param ... `r roxy_dots()`
#'
#' @return a [data.table][data.table::data.table-package] containing
#'   partial dependence values for the specified variable(s)
#'   and, if relevant, at the specified prediction horizon(s).
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
                        pred_type = NULL,
                        expand_grid = TRUE,
                        prob_values = c(0.025, 0.50, 0.975),
                        prob_labels = c('lwr', 'medn', 'upr'),
                        boundary_checks = TRUE,
                        n_thread = 0,
                        verbose_progress = FALSE,
                        ...){

 check_dots(list(...), orsf_pd_oob)

 orsf_dependence(object = object,
                 pred_spec = pred_spec,
                 pd_data = NULL,
                 pred_horizon = pred_horizon,
                 pred_type = pred_type,
                 expand_grid = expand_grid,
                 prob_values = prob_values,
                 prob_labels = prob_labels,
                 boundary_checks = boundary_checks,
                 n_thread = n_thread,
                 verbose_progress = verbose_progress,
                 oobag = TRUE,
                 type_output = 'smry')

}

#' @rdname orsf_pd_oob
#' @export
orsf_pd_inb <- function(object,
                        pred_spec,
                        pred_horizon = NULL,
                        pred_type = NULL,
                        expand_grid = TRUE,
                        prob_values = c(0.025, 0.50, 0.975),
                        prob_labels = c('lwr', 'medn', 'upr'),
                        boundary_checks = TRUE,
                        n_thread = 0,
                        verbose_progress = FALSE,
                        ...){

 check_dots(list(...), orsf_pd_inb)

 if(is.null(object$data))
  stop("no data were found in object. ",
       "did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)

 orsf_dependence(object = object,
                 pred_spec = pred_spec,
                 pd_data = object$data,
                 pred_horizon = pred_horizon,
                 pred_type = pred_type,
                 expand_grid = expand_grid,
                 prob_values = prob_values,
                 prob_labels = prob_labels,
                 boundary_checks = boundary_checks,
                 n_thread = n_thread,
                 verbose_progress = verbose_progress,
                 oobag = FALSE,
                 type_output = 'smry')

}

#' @rdname orsf_pd_oob
#' @export
orsf_pd_new <- function(object,
                        pred_spec,
                        new_data,
                        pred_horizon = NULL,
                        pred_type = NULL,
                        na_action = 'fail',
                        expand_grid = TRUE,
                        prob_values = c(0.025, 0.50, 0.975),
                        prob_labels = c('lwr', 'medn', 'upr'),
                        boundary_checks = TRUE,
                        n_thread = 0,
                        verbose_progress = FALSE,
                        ...){

 check_dots(list(...), orsf_pd_new)

 orsf_dependence(object = object,
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
                 verbose_progress = verbose_progress,
                 oobag = FALSE,
                 type_output = 'smry')

}


#' Individual Conditional Expectations
#'
#' Compute individual conditional expectations for an
#'  oblique random forest. `r roxy_ice_explain()`
#'  `r roxy_pd_oob_explain('individual conditional expectations')`
#'
#' @inheritParams orsf_pd_oob
#'
#' @return a [data.table][data.table::data.table-package] containing
#'   individual conditional expectations for the specified variable(s)
#'   and, if relevant, at the specified prediction horizon(s).
#'
#' @includeRmd Rmd/orsf_ice_examples.Rmd
#'
#' @export
#'
#'
orsf_ice_oob <- function(object,
                         pred_spec,
                         pred_horizon = NULL,
                         pred_type = NULL,
                         expand_grid = TRUE,
                         boundary_checks = TRUE,
                         n_thread = 0,
                         verbose_progress = FALSE,
                         ...){

 check_dots(list(...), orsf_ice_oob)

 orsf_dependence(object = object,
                 pred_spec = pred_spec,
                 pd_data = NULL,
                 pred_horizon = pred_horizon,
                 pred_type = pred_type,
                 expand_grid = expand_grid,
                 boundary_checks = boundary_checks,
                 n_thread = n_thread,
                 verbose_progress = verbose_progress,
                 oobag = TRUE,
                 type_output = 'ice')

}

#' @rdname orsf_ice_oob
#' @export
orsf_ice_inb <- function(object,
                         pred_spec,
                         pred_horizon = NULL,
                         pred_type = NULL,
                         expand_grid = TRUE,
                         boundary_checks = TRUE,
                         n_thread = 0,
                         verbose_progress = FALSE,
                         ...){

 check_dots(list(...), orsf_ice_oob)

 if(is.null(object$data))
  stop("no data were found in object. ",
       "did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)

 orsf_dependence(object = object,
                 pred_spec = pred_spec,
                 pd_data = object$data,
                 pred_horizon = pred_horizon,
                 pred_type = pred_type,
                 expand_grid = expand_grid,
                 boundary_checks = boundary_checks,
                 n_thread = n_thread,
                 verbose_progress = verbose_progress,
                 oobag = FALSE,
                 type_output = 'ice')

}

#' @rdname orsf_ice_oob
#' @export
orsf_ice_new <- function(object,
                         pred_spec,
                         new_data,
                         pred_horizon = NULL,
                         pred_type = NULL,
                         na_action = 'fail',
                         expand_grid = TRUE,
                         boundary_checks = TRUE,
                         n_thread = 0,
                         verbose_progress = FALSE,
                         ...){

 check_dots(list(...), orsf_ice_new)

 orsf_dependence(object = object,
                 pred_spec = pred_spec,
                 pd_data = new_data,
                 pred_horizon = pred_horizon,
                 pred_type = pred_type,
                 na_action = na_action,
                 expand_grid = expand_grid,
                 boundary_checks = boundary_checks,
                 n_thread = n_thread,
                 verbose_progress = verbose_progress,
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

orsf_dependence <- function(object,
                            pd_data,
                            pred_spec,
                            pred_horizon,
                            pred_type,
                            na_action = NULL,
                            expand_grid,
                            prob_values = NULL,
                            prob_labels = NULL,
                            boundary_checks,
                            n_thread,
                            verbose_progress,
                            oobag,
                            type_output){

 check_arg_is(object, arg_name = 'object', expected_class = 'ObliqueForest')

 if(oobag && is.null(object$data))
  stop("no data were found in object. ",
       "did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)

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
                           verbose_progress = verbose_progress,
                           oobag = oobag,
                           type_output = type_output)


}
