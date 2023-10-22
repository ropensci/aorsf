
#' null operator (copied from rlang)

`%||%` <-  function (x, y) {
 if (is.null(x))
  y
 else x
}

#' helper for guessing pred_horizon input
#'
#' @param object 'orsf_fit' object
#' @param pred_horizon NULL or a user's specified pred_horizon
#'
#' @return
#'  - if the user gave a pred_horizon, return that.
#'  - else if the object has a pred_horizon, return that
#'  - else throw an error
#'
#' @noRd

infer_pred_horizon <- function(object, pred_type, pred_horizon){

 check_arg_is(object, 'object', 'orsf_fit')

 if(pred_type %in% c("mort", "leaf")){
  # value of pred_horizon does not matter for these types of prediction
  pred_horizon <- 1
 }

 # see if it was previously specified
 if(is.null(pred_horizon)) pred_horizon <- object$pred_horizon

 # throw error if pred_type requires pred_horizon
 if(is.null(pred_horizon)){
  stop("pred_horizon was not specified and could not be found in object.",
       call. = FALSE)
 }


 pred_horizon

}


#' helper for guessing outcome type
#'
#' @param names_y_data character vector of outcome names
#' @param data dataset containing outcomes
#'
#' @return character value: 'survival', 'regression' or 'classification'
#'
#' @examples
#'
#' infer_outcome_type('bili', pbc_orsf)
#' infer_outcome_type('sex', pbc_orsf)
#' infer_outcome_type(c('time', 'status'), pbc_orsf)
#' infer_outcome_type(Surv(pbc_orsf$time, pbc_orsf$status), pbc_orsf)
#'
#' @noRd
infer_outcome_type <- function(names_y_data, data){

 if(length(names_y_data) > 2){
  stop("formula should have at most two variables as the response",
       call. = FALSE)
 }

 if(length(names_y_data) == 2) {
  return("survival")
 }

 if(is.factor(data[[names_y_data]])){
  return("classification")
 } else if(inherits(data[[names_y_data]], 'Surv')) {
  return("survival")
 } else {
  return("regression")
 }

 stop("could not infer outcome type", call. = FALSE)

}


infer_orsf_args <- function(x,
                            y,
                            w = rep(1, nrow(x)),
                            ...,
                            object = NULL){

 .dots <- list(...)

 control <- .dots$control %||%
  get_control(object) %||%
  orsf_control_fast()

 n_tree = .dots$n_tree %||%
  get_n_tree(object) %||%
  500L

 tree_type = .dots$tree_type %||%
  get_tree_type(object) %||%
  'survival'

 split_rule <- .dots$split_rule %||%
  get_split_rule(object) %||%
  'logrank'

 split_min_stat <- .dots$split_min_stat %||%
  get_split_min_stat(object) %||%
  switch(split_rule, "logrank" = 3.841459, "cstat" = 0.50)

 mtry <- .dots$mtry %||%
  get_mtry(object) %||%
  ceiling(sqrt(ncol(x)))

 oobag_pred_type <- get_oobag_pred_type(object) %||% "surv"
 oobag_pred <- oobag_pred_type != 'none'


 pred_horizon <- get_oobag_pred_horizon(object) %||%
  if(tree_type == 'survival') stats::median(y[, 1]) else 1

 type_oobag_eval <- get_type_oobag_eval(object) %||%
  if(oobag_pred) 'cstat' else 'none'

 vi_type <- .dots$vi_type %||% get_importance(object) %||% "anova"

 pd_type <- .dots$pd_type %||% 'none'

 list(
  x = x,
  y = y,
  w = w,
  tree_type_R = switch(tree_type,
                       'classification' = 1,
                       'regression'= 2,
                       'survival' = 3),
  tree_seeds = .dots$tree_seeds %||% get_tree_seeds(object) %||% 329,
  loaded_forest = object$forest %||% list(),
  n_tree = n_tree,
  mtry = mtry,
  sample_with_replacement = .dots$sample_with_replacement %||%
   get_sample_with_replacement(object) %||%
   TRUE,
  sample_fraction = .dots$sample_fraction %||%
   get_sample_fraction(object) %||%
   0.632,
  vi_type_R = switch(vi_type,
                     "none" = 0,
                     "negate" = 1,
                     "permute" = 2,
                     "anova" = 3),
  vi_max_pvalue = .dots$vi_max_pvalue %||%
   get_vi_max_pvalue(object) %||%
   0.01,
  leaf_min_events = .dots$leaf_min_events %||%
   get_leaf_min_events(object) %||%
   1,
  leaf_min_obs = .dots$leaf_min_obs %||%
   get_leaf_min_obs(object) %||%
   5,
  split_rule_R = switch(split_rule, "logrank" = 1, "cstat" = 2),
  split_min_events = .dots$split_min_event %||%
   get_split_min_events(object) %||%
   5,
  split_min_obs = .dots$split_min_obs %||%
   get_split_min_obs(object) %||%
   10,
  split_min_stat = .dots$split_min_stat %||%
   get_split_min_stat(object) %||%
   NA_real_,
  split_max_cuts = .dots$split_max_cuts %||%
   get_n_split(object) %||%
   5,
  split_max_retry = .dots$split_max_retry %||%
   get_n_retry(object) %||%
   3,
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
  lincomb_df_target = control$lincomb_df_target %||% mtry,
  lincomb_ties_method = switch(tolower(control$lincomb_ties_method),
                               'breslow' = 0,
                               'efron'   = 1),
  pred_type_R = switch(oobag_pred_type,
                       "none" = 0,
                       "risk" = 1,
                       "surv" = 2,
                       "chf"  = 3,
                       "mort" = 4,
                       "leaf" = 8),
  pred_mode = .dots$pred_mode %||% FALSE,
  pred_aggregate = .dots$pred_aggregate %||% oobag_pred_type != 'leaf',
  pred_horizon = pred_horizon,
  oobag = oobag_pred,
  oobag_R_function = .dots$oobag_R_function %||%
   get_f_oobag_eval(object) %||%
   function(x) x,
  oobag_eval_type_R = switch(type_oobag_eval,
                             'none' = 0,
                             'cstat' = 1,
                             'user' = 2),
  oobag_eval_every = .dots$oobag_eval_every %||%
   get_oobag_eval_every(object) %||%
   n_tree,
  pd_type_R = switch(pd_type, "none" = 0L, "smry" = 1L, "ice" = 2L),
  pd_x_vals = .dots$pd_x_vals %||% list(matrix(0, ncol=0, nrow=0)),
  pd_x_cols = .dots$pd_x_cols %||% list(matrix(0, ncol=0, nrow=0)),
  pd_probs  = .dots$pd_probs %||% 0,
  n_thread  = .dots$n_thread %||% get_n_thread(object) %||% 1,
  write_forest = .dots$write_forest %||% TRUE,
  run_forest   = .dots$run_forest %||% TRUE,
  verbosity    = .dots$verbosity %||%
   get_verbose_progress(object) %||%
   FALSE
 )

}
