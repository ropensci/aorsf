

#' Compute predictions using ORSF
#'
#' Predicted risk, survival, hazard, or mortality from an ORSF model.
#'
#' @param object (*orsf_fit*) a trained oblique random survival forest
#'   (see [orsf]).
#'
#' @param new_data a `r roxy_data_allowed()` to compute predictions in.
#'
#' @param pred_horizon (_double_) a value or vector indicating the time(s)
#'   that predictions will be calibrated to. E.g., if you were predicting
#'   risk of incident heart failure within the next 10 years, then
#'   `pred_horizon = 10`. `pred_horizon` can be `NULL` if `pred_type` is
#'   `'mort'`, since mortality predictions are aggregated over all
#'   event times
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
#'   - `r roxy_na_action_pass("new_data")`
#'   - `r roxy_na_action_omit("new_data")`
#'   - `r roxy_na_action_impute_meanmode('new_data')`. To clarify,
#'     the mean and mode used to impute missing values are from the
#'     training data of `object`, not from `new_data`.
#'
#' @param boundary_checks (_logical_) if `TRUE`, `pred_horizon` will be
#'  checked to make sure the requested values are less than the maximum
#'  observed time in `object`'s training data. If `FALSE`, these checks
#'  are skipped.
#'
#' @param n_thread `r roxy_n_thread_header("computing predictions")`
#'
#' @param pred_aggregate (_logical_) If `TRUE` (the default), predictions
#'   will be aggregated over all trees by taking the mean. If `FALSE`, the
#'   returned output will contain one row per observation and one column
#'   for each tree. If the length of `pred_horizon` is two or more and
#'   `pred_aggregate` is `FALSE`, then the result will be a list of such
#'   matrices, with the i'th item in the list corresponding to the i'th
#'   value of `pred_horizon`.
#'
#' @inheritParams orsf
#'
#' @param ... `r roxy_dots()`
#'
#' @return a `matrix` of predictions. Column `j` of the matrix corresponds
#'   to value `j` in `pred_horizon`. Row `i` of the matrix corresponds to
#'   row `i` in `new_data`.
#'
#' @details
#'
#' `new_data` must have the same columns with equivalent types as the data
#'   used to train `object`. Also, factors in `new_data` must not have levels
#'   that were not in the data used to train `object`.
#'
#' `pred_horizon` values should not exceed the maximum follow-up time in
#'   `object`'s training data, but if you truly want to do this, set
#'   `boundary_checks = FALSE` and you can use a `pred_horizon` as large
#'   as you want. Note that predictions beyond the maximum follow-up time
#'   in the `object`'s training data are equal to predictions at the
#'   maximum follow-up time, because `aorsf` does not estimate survival
#'   beyond its maximum observed time.
#'
#' If unspecified, `pred_horizon` may be automatically specified as the value
#'   used for `oobag_pred_horizon` when `object` was created (see [orsf]).
#'
#'
#' @export
#'
#' @includeRmd Rmd/orsf_predict_examples.Rmd
#'
predict.orsf_fit <- function(object,
                             new_data,
                             pred_horizon = NULL,
                             pred_type = 'risk',
                             na_action = 'fail',
                             boundary_checks = TRUE,
                             n_thread = 1,
                             verbose_progress = FALSE,
                             pred_aggregate = TRUE,
                             ...){

 # catch any arguments that didn't match and got relegated to ...
 # these arguments are mistaken input names since ... isn't used.
 check_dots(list(...), .f = predict.orsf_fit)

 names_x_data <- intersect(get_names_x(object), names(new_data))

 if(pred_type %in% c('leaf', 'mort') && !is.null(pred_horizon)){

  extra_text <- if(length(pred_horizon)>1){
   " Predictions at each value of pred_horizon will be identical."
  } else {
   ""
  }

  warning("pred_horizon does not impact predictions",
          " when pred_type is '", pred_type, "'.",
          extra_text, call. = FALSE)

 }

 pred_horizon <- infer_pred_horizon(object, pred_type, pred_horizon)

 check_predict(object = object,
               new_data = new_data,
               pred_horizon = pred_horizon,
               pred_type = pred_type,
               na_action = na_action,
               boundary_checks = boundary_checks)

 if(length(pred_horizon) > 1 && !pred_aggregate){

  results <- lapply(
   X = pred_horizon,
   FUN = function(t){
    predict.orsf_fit(object = object,
                     new_data = new_data,
                     pred_horizon = t,
                     pred_type = pred_type,
                     na_action = na_action,
                     boundary_checks = boundary_checks,
                     n_thread = n_thread,
                     verbose_progress = verbose_progress,
                     pred_aggregate = pred_aggregate)
   }
  )

  return(simplify2array(results))

 }

 pred_horizon_order <- order(pred_horizon)
 pred_horizon_ordered <- pred_horizon[pred_horizon_order]

 cc <- which(
  stats::complete.cases(select_cols(new_data, names_x_data))
 )

 check_complete_cases(cc, na_action, nrow(new_data))

 if(na_action == 'impute_meanmode'){

  new_data <- data_impute(new_data,
                          cols = get_names_x(object),
                          values = c(as.list(get_means(object)),
                                     as.list(get_modes(object))))

  cc <- collapse::seq_row(new_data)

 }

 if(is.null(pred_horizon) && !(pred_type %in% c('mort', 'leaf'))){
  stop("pred_horizon must be specified for ",
       pred_type, " predictions.", call. = FALSE)
 }

 x_new <- prep_x_from_orsf(object, data = new_data[cc, ])

 # control <- get_control(object)

 args <- infer_orsf_args(x = x_new,
                         vi_type = 'none',
                         object = object,
                         pred_type = pred_type,
                         pred_aggregate = pred_aggregate,
                         pred_horizon = pred_horizon_ordered,
                         oobag_pred = FALSE,
                         pred_mode = TRUE,
                         write_forest = FALSE,
                         run_forest = TRUE)

 orsf_out <- do.call(orsf_cpp, args = args)

 out_values <- orsf_out$pred_new

 if(na_action == "pass"){

  out <- matrix(nrow = nrow(new_data),
                ncol = ncol(out_values))

  out[cc, ] <- out_values

 } else {

  out <- out_values

 }

 if(pred_type == "leaf" || !pred_aggregate) return(out)

 if(get_tree_type(object) == 'survival'){
  # output in the same order as pred_horizon
  out <- out[, order(pred_horizon_order), drop = FALSE]
 }

 out

}












# old code in case the infer function fails me:
# args_tmp <- list(x = x_new,
#                  y = matrix(1, ncol=2),
#                  w = rep(1, nrow(x_new)),
#                  tree_type_R = get_tree_type(object),
#                  tree_seeds = get_tree_seeds(object),
#                  loaded_forest = object$forest,
#                  n_tree = get_n_tree(object),
#                  mtry = get_mtry(object),
#                  sample_with_replacement = get_sample_with_replacement(object),
#                  sample_fraction = get_sample_fraction(object),
#                  vi_type_R = 0,
#                  vi_max_pvalue = get_importance_max_pvalue(object),
#                  oobag_R_function = get_f_oobag_eval(object),
#                  leaf_min_events = get_leaf_min_events(object),
#                  leaf_min_obs = get_leaf_min_obs(object),
#                  split_rule_R = switch(get_split_rule(object),
#                                        "logrank" = 1,
#                                        "cstat" = 2),
#                  split_min_events = get_split_min_events(object),
#                  split_min_obs = get_split_min_obs(object),
#                  split_min_stat = get_split_min_stat(object),
#                  split_max_cuts = get_n_split(object),
#                  split_max_retry = get_n_retry(object),
#                  lincomb_R_function = control$lincomb_R_function,
#                  lincomb_type_R = switch(control$lincomb_type,
#                                          'glm' = 1,
#                                          'random' = 2,
#                                          'net' = 3,
#                                          'custom' = 4),
#                  lincomb_eps = control$lincomb_eps,
#                  lincomb_iter_max = control$lincomb_iter_max,
#                  lincomb_scale = control$lincomb_scale,
#                  lincomb_alpha = control$lincomb_alpha,
#                  lincomb_df_target = control$lincomb_df_target,
#                  lincomb_ties_method = switch(tolower(control$lincomb_ties_method),
#                                               'breslow' = 0,
#                                               'efron'   = 1),
#                  pred_type_R = switch(pred_type,
#                                       "risk" = 1,
#                                       "surv" = 2,
#                                       "chf"  = 3,
#                                       "mort" = 4,
#                                       "leaf" = 8),
#                  pred_mode = TRUE,
#                  pred_aggregate = pred_aggregate,
#                  pred_horizon = pred_horizon_ordered,
#                  oobag = FALSE,
#                  oobag_eval_type_R = 0,
#                  oobag_eval_every = get_n_tree(object),
#                  pd_type_R = 0,
#                  pd_x_vals = list(matrix(0, ncol=1, nrow=1)),
#                  pd_x_cols = list(matrix(1L, ncol=1, nrow=1)),
#                  pd_probs = c(0),
#                  n_thread = n_thread,
#                  write_forest = FALSE,
#                  run_forest = TRUE,
#                  verbosity = as.integer(verbose_progress))
#
# checkout <- c()
#
# for(i in names(args)){
#  print(i)
#  if(!is.list(args[[i]]) && !is.function(args[[i]])){
#   if(!all(args[[i]] == args_tmp[[i]]))
#    checkout <- c(checkout, i)
#  }
# }
# browser()
# orsf_out <- orsf_cpp(x = x_new,
#                      y = matrix(1, ncol=2),
#                      w = rep(1, nrow(x_new)),
#                      tree_type_R = get_tree_type(object),
#                      tree_seeds = get_tree_seeds(object),
#                      loaded_forest = object$forest,
#                      n_tree = get_n_tree(object),
#                      mtry = get_mtry(object),
#                      sample_with_replacement = get_sample_with_replacement(object),
#                      sample_fraction = get_sample_fraction(object),
#                      vi_type_R = 0,
#                      vi_max_pvalue = get_importance_max_pvalue(object),
#                      oobag_R_function = get_f_oobag_eval(object),
#                      leaf_min_events = get_leaf_min_events(object),
#                      leaf_min_obs = get_leaf_min_obs(object),
#                      split_rule_R = switch(get_split_rule(object),
#                                            "logrank" = 1,
#                                            "cstat" = 2),
#                      split_min_events = get_split_min_events(object),
#                      split_min_obs = get_split_min_obs(object),
#                      split_min_stat = get_split_min_stat(object),
#                      split_max_cuts = get_n_split(object),
#                      split_max_retry = get_n_retry(object),
#                      lincomb_R_function = control$lincomb_R_function,
#                      lincomb_type_R = switch(control$lincomb_type,
#                                              'glm' = 1,
#                                              'random' = 2,
#                                              'net' = 3,
#                                              'custom' = 4),
#                      lincomb_eps = control$lincomb_eps,
#                      lincomb_iter_max = control$lincomb_iter_max,
#                      lincomb_scale = control$lincomb_scale,
#                      lincomb_alpha = control$lincomb_alpha,
#                      lincomb_df_target = control$lincomb_df_target,
#                      lincomb_ties_method = switch(tolower(control$lincomb_ties_method),
#                                                   'breslow' = 0,
#                                                   'efron'   = 1),
#                      pred_type_R = switch(pred_type,
#                                           "risk" = 1,
#                                           "surv" = 2,
#                                           "chf"  = 3,
#                                           "mort" = 4,
#                                           "leaf" = 8),
#                      pred_mode = TRUE,
#                      pred_aggregate = pred_aggregate,
#                      pred_horizon = pred_horizon_ordered,
#                      oobag = FALSE,
#                      oobag_eval_type_R = 0,
#                      oobag_eval_every = get_n_tree(object),
#                      pd_type_R = 0,
#                      pd_x_vals = list(matrix(0, ncol=1, nrow=1)),
#                      pd_x_cols = list(matrix(1L, ncol=1, nrow=1)),
#                      pd_probs = c(0),
#                      n_thread = n_thread,
#                      write_forest = FALSE,
#                      run_forest = TRUE,
#                      verbosity = as.integer(verbose_progress))

