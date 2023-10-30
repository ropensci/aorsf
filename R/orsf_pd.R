


#' ORSF partial dependence
#'
#' Compute partial dependence for an ORSF model.
#' `r roxy_pd_explain()`
#' `r roxy_pd_oob_explain('partial dependence')`
#'
#' @inheritParams predict.orsf_fit
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

 pred_horizon <- infer_pred_horizon(object, pred_type, pred_horizon)

 # make a visible binding for CRAN
 id_variable = NULL

 if(is.null(prob_values)) prob_values <- c(0.025, 0.50, 0.975)
 if(is.null(prob_labels)) prob_labels <- c('lwr', 'medn', 'upr')

 check_pd_inputs(object          = object,
                 pred_spec       = pred_spec,
                 expand_grid     = expand_grid,
                 prob_values     = prob_values,
                 prob_labels     = prob_labels,
                 oobag           = oobag,
                 boundary_checks = boundary_checks,
                 new_data        = pd_data,
                 pred_type       = pred_type,
                 pred_horizon    = pred_horizon,
                 na_action       = na_action)

 if(oobag && is.null(object$data))
  stop("no data were found in object. ",
       "did you use attach_data = FALSE when ",
       "running orsf()?", call. = FALSE)

 if(oobag) pd_data <- object$data

 if(is.null(pred_horizon) && pred_type != 'mort'){
  stop("pred_horizon must be specified for ",
       pred_type, " predictions.", call. = FALSE)
 }

 names_x_data <- intersect(get_names_x(object), names(pd_data))

 cc <- which(stats::complete.cases(select_cols(pd_data, names_x_data)))

 if(na_action == 'pass'){
  stop("na_action = 'pass' is not supported for individual conditional ",
       "expectation. Please use na_action = 'omit' or na_action = 'fail'.",
       call. = FALSE)
 }

 check_complete_cases(cc, na_action, nrow(pd_data))

 x_new <- prep_x_from_orsf(object, data = pd_data[cc, ])

 # the values in pred_spec need to be centered & scaled to match x_new,
 # which is also centered and scaled
 means <- get_means(object)
 standard_deviations <- get_standard_deviations(object)

 for(i in intersect(names(means), names(pred_spec))){
  pred_spec[[i]] <- (pred_spec[[i]] - means[i]) / standard_deviations[i]
 }

 pred_type_R <- switch(pred_type,
                       "risk" = 1,
                       "surv" = 2,
                       "chf"  = 3,
                       "mort" = 4)

 fi <- get_fctr_info(object)

 if(expand_grid){

  if(!is.data.frame(pred_spec))
   pred_spec <- expand.grid(pred_spec, stringsAsFactors = TRUE)

  for(i in seq_along(fi$cols)){

   ii <- fi$cols[i]

   if(is.character(pred_spec[[ii]]) && !fi$ordr[i]){

    pred_spec[[ii]] <- factor(pred_spec[[ii]], levels = fi$lvls[[ii]])

   }

  }

  check_new_data_fctrs(new_data  = pred_spec,
                       names_x   = get_names_x(object),
                       fi_ref    = fi,
                       label_new = "pred_spec")

  pred_spec_new <- ref_code(x_data = pred_spec,
                            fi = get_fctr_info(object),
                            names_x_data = names(pred_spec))

  x_cols <- list(match(names(pred_spec_new), colnames(x_new)) - 1)

  pred_spec_new <- list(as.matrix(pred_spec_new))

  pd_bind <- list(pred_spec)

 } else {

  pred_spec_new <- pd_bind <- x_cols <- list()

  for(i in seq_along(pred_spec)){

   pred_spec_new[[i]]  <- as.data.frame(pred_spec[i])
   pd_name <- names(pred_spec)[i]

   pd_bind[[i]] <- data.frame(
    variable = pd_name,
    value = rep(NA_real_, length(pred_spec[[i]])),
    level = rep(NA_character_, length(pred_spec[[i]]))
   )

   if(pd_name %in% fi$cols) {

    pd_bind[[i]]$level <- as.character(pred_spec[[i]])

    pred_spec_new[[i]] <- ref_code(pred_spec_new[[i]],
                                   fi = fi,
                                   names_x_data = pd_name)

   } else {

    pd_bind[[i]]$value <- pred_spec[[i]]

   }

   x_cols[[i]] <- match(names(pred_spec_new[[i]]), colnames(x_new)) - 1
   pred_spec_new[[i]] <- as.matrix(pred_spec_new[[i]])

  }

 }

 control <- get_control(object)

 pred_horizon_order <- order(pred_horizon)
 pred_horizon_ordered <- pred_horizon[pred_horizon_order]

 # results <- list()
 #
 # for(i in seq_along(pred_spec_new)){
 #
 #  results_i <- list()
 #
 #  x_pd <- x_new
 #
 #  for(j in seq(nrow(pred_spec_new[[i]]))){
 #
 #   x_pd[, x_cols[[i]]] <- pred_spec_new[[i]][j, ]
 #
 #   results_i[[j]] <- orsf_cpp(
 #    x = x_pd,
 #    y = matrix(1, ncol=2),
 #    w = rep(1, nrow(x_new)),
 #    tree_type_R = get_tree_type(object),
 #    tree_seeds = get_tree_seeds(object),
 #    loaded_forest = object$forest,
 #    n_tree = get_n_tree(object),
 #    mtry = get_mtry(object),
 #    sample_with_replacement = get_sample_with_replacement(object),
 #    sample_fraction = get_sample_fraction(object),
 #    vi_type_R = 0,
 #    vi_max_pvalue = get_vi_max_pvalue(object),
 #    oobag_R_function = get_f_oobag_eval(object),
 #    leaf_min_events = get_leaf_min_events(object),
 #    leaf_min_obs = get_leaf_min_obs(object),
 #    split_rule_R = switch(get_split_rule(object),
 #                          "logrank" = 1,
 #                          "cstat" = 2),
 #    split_min_events = get_split_min_events(object),
 #    split_min_obs = get_split_min_obs(object),
 #    split_min_stat = get_split_min_stat(object),
 #    split_max_cuts = get_n_split(object),
 #    split_max_retry = get_n_retry(object),
 #    lincomb_R_function = control$lincomb_R_function,
 #    lincomb_type_R = switch(control$lincomb_type,
 #                            'glm' = 1,
 #                            'random' = 2,
 #                            'net' = 3,
 #                            'custom' = 4),
 #    lincomb_eps = control$lincomb_eps,
 #    lincomb_iter_max = control$lincomb_iter_max,
 #    lincomb_scale = control$lincomb_scale,
 #    lincomb_alpha = control$lincomb_alpha,
 #    lincomb_df_target = control$lincomb_df_target,
 #    lincomb_ties_method = switch(tolower(control$lincomb_ties_method),
 #                                 'breslow' = 0,
 #                                 'efron'   = 1),
 #    pred_type_R = pred_type_R,
 #    pred_mode = TRUE,
 #    pred_aggregate = TRUE,
 #    pred_horizon = pred_horizon_ordered,
 #    oobag = oobag,
 #    oobag_eval_type_R = 0,
 #    oobag_eval_every = get_n_tree(object),
 #    pd_type_R = 0,
 #    pd_x_vals = list(matrix(0, ncol=1, nrow=1)),
 #    pd_x_cols = list(matrix(1L, ncol=1, nrow=1)),
 #    pd_probs = c(0),
 #    n_thread = n_thread,
 #    write_forest = FALSE,
 #    run_forest = TRUE,
 #    verbosity = 0)$pred_new
 #
 #  }
 #
 #  if(type_output == 'smry'){
 #   results_i <- lapply(
 #    results_i,
 #    function(x) {
 #     apply(x, 2, function(x_col){
 #      as.numeric(
 #       c(mean(x_col, na.rm = TRUE),
 #         quantile(x_col, probs = prob_values, na.rm = TRUE))
 #      )
 #     })
 #    }
 #   )
 #  }
 #
 #
 #  results[[i]] <- results_i
 #
 # }
 #
 # pd_vals <- results

 orsf_out <- orsf_cpp(x = x_new,
                      y = matrix(1, ncol=2),
                      w = rep(1, nrow(x_new)),
                      tree_type_R = get_tree_type(object),
                      tree_seeds = get_tree_seeds(object),
                      loaded_forest = object$forest,
                      n_tree = get_n_tree(object),
                      mtry = get_mtry(object),
                      sample_with_replacement = get_sample_with_replacement(object),
                      sample_fraction = get_sample_fraction(object),
                      vi_type_R = 0,
                      vi_max_pvalue = get_vi_max_pvalue(object),
                      oobag_R_function = get_f_oobag_eval(object),
                      leaf_min_events = get_leaf_min_events(object),
                      leaf_min_obs = get_leaf_min_obs(object),
                      split_rule_R = switch(get_split_rule(object),
                                            "logrank" = 1,
                                            "cstat" = 2),
                      split_min_events = get_split_min_events(object),
                      split_min_obs = get_split_min_obs(object),
                      split_min_stat = get_split_min_stat(object),
                      split_max_cuts = get_n_split(object),
                      split_max_retry = get_n_retry(object),
                      lincomb_R_function = control$lincomb_R_function,
                      lincomb_type_R = switch(control$lincomb_type,
                                              'glm' = 1,
                                              'random' = 2,
                                              'net' = 3,
                                              'custom' = 4),
                      lincomb_eps = control$lincomb_eps,
                      lincomb_iter_max = control$lincomb_iter_max,
                      lincomb_scale = control$lincomb_scale,
                      lincomb_alpha = control$lincomb_alpha,
                      lincomb_df_target = control$lincomb_df_target,
                      lincomb_ties_method = switch(tolower(control$lincomb_ties_method),
                                                   'breslow' = 0,
                                                   'efron'   = 1),
                      pred_type_R = pred_type_R,
                      pred_mode = FALSE,
                      pred_aggregate = TRUE,
                      pred_horizon = pred_horizon,
                      oobag = oobag,
                      oobag_eval_type_R = 0,
                      oobag_eval_every = get_n_tree(object),
                      pd_type_R = switch(type_output,
                                         "smry" = 1L,
                                         "ice" = 2L),
                      pd_x_vals = pred_spec_new,
                      pd_x_cols = x_cols,
                      pd_probs = prob_values,
                      n_thread = n_thread,
                      write_forest = FALSE,
                      run_forest = TRUE,
                      verbosity = 0)

 pd_vals <- orsf_out$pd_values

 for(i in seq_along(pd_vals)){

  pd_bind[[i]]$id_variable <- seq(nrow(pd_bind[[i]]))

  for(j in seq_along(pd_vals[[i]])){

   pd_vals[[i]][[j]] <- matrix(pd_vals[[i]][[j]],
                               nrow=length(pred_horizon),
                               byrow = T)

   rownames(pd_vals[[i]][[j]]) <- pred_horizon

   if(type_output=='smry')
    colnames(pd_vals[[i]][[j]]) <- c('mean', prob_labels)
   else
    colnames(pd_vals[[i]][[j]]) <- c(paste(1:nrow(x_new)))

   ph <- rownames(pd_vals[[i]][[j]])

   pd_vals[[i]][[j]] <- as.data.frame(pd_vals[[i]][[j]])

   rownames(pd_vals[[i]][[j]]) <- NULL

   pd_vals[[i]][[j]][['pred_horizon']] <- ph

   if(type_output == 'ice'){

    pd_vals[[i]][[j]] <- melt_aorsf(
     data = pd_vals[[i]][[j]],
     id.vars = 'pred_horizon',
     variable.name = 'id_row',
     value.name = 'pred',
     measure.vars = setdiff(names(pd_vals[[i]][[j]]), 'pred_horizon'))

   }


  }

  pd_vals[[i]] <- rbindlist(pd_vals[[i]], idcol = 'id_variable')

  # this seems awkward but the reason I convert back to data.frame
  # here is to avoid a potential memory leak from forder & bmerge.
  # I have no idea why this memory leak may be occurring but it does
  # not if I apply merge.data.frame instead of merge.data.table
  pd_vals[[i]] <- merge(as.data.frame(pd_vals[[i]]),
                        as.data.frame(pd_bind[[i]]),
                        by = 'id_variable')

 }

 out <- rbindlist(pd_vals)

 # # missings may occur when oobag=TRUE and n_tree is small
 # if(type_output == 'ice') {
 #  out <- collapse::na_omit(out, cols = 'pred')
 # }

 ids <- c('id_variable')

 if(type_output == 'ice') ids <- c(ids, 'id_row')

 mid <- setdiff(names(out), c(ids, 'mean', prob_labels, 'pred'))

 end <- setdiff(names(out), c(ids, mid))

 setcolorder(out, neworder = c(ids, mid, end))

 out[, pred_horizon := as.numeric(pred_horizon)]

 # not needed for summary
 if(type_output == 'smry')
  out[, id_variable := NULL]

 # not needed for mort
 if(pred_type == 'mort')
  out[, pred_horizon := NULL]

 # put data back into original scale
 for(j in intersect(names(means), names(pred_spec))){

  if(j %in% names(out)){

   var_index <- collapse::seq_row(out)
   var_value <- (out[[j]] * standard_deviations[j]) + means[j]
   var_name  <- j

  } else {

   var_index <- out$variable %==% j
   var_value <- (out$value[var_index] * standard_deviations[j]) + means[j]
   var_name  <- 'value'

  }

  set(out, i = var_index, j = var_name, value = var_value)

 }

 # silent print after modify in place
 out[]

 out

}


pd_list_split <- function(x_vals, x_cols){

 x_vals_out <- x_cols_out <- vector(mode = 'list')
 counter <- 1

 for(i in seq_along(x_vals)){

  x_vals_split <- split(x_vals[[i]], row(x_vals[[i]]))

  for(j in seq_along(x_vals_split)){

   x_vals_out[[counter]] <- matrix(x_vals_split[[j]],
                                   ncol = ncol(x_vals[[i]]),
                                   nrow = 1)
   colnames(x_vals_out[[counter]]) <- colnames(x_vals[[i]])

   x_cols_out[[counter]] <- x_cols[[i]]

   counter <- counter + 1

  }

 }

 list(
  x_vals = x_vals_out,
  x_cols = x_cols_out
 )

}

