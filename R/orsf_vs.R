
#' Variable selection
#'
#' @inheritParams predict.orsf_fit
#' @param n_predictor_min (*integer*) the minimum number of predictors allowed
#' @param verbose_progress (*logical*) not implemented yet. Should progress be printed to the console?
#'
#' @return a [data.table][data.table::data.table-package] with four columns:
#'   - *n_predictors*: the number of predictors used
#'   - *stat_value*: the out-of-bag statistic
#'   - *predictors_included*: the names of the predictors included
#'   - *predictor_dropped*: the predictor selected to be dropped
#'
#' @details
#'
#' `tree_seeds` should be specified in `object` so that each successive run
#'   of `orsf` will be evaluated in the same out-of-bag samples as the initial
#'   run.
#'
#' @export
#'
#' @examples
#'
#' object <- orsf(formula = time + status ~ .,
#'                data = pbc_orsf,
#'                n_tree = 25,
#'                importance = 'anova',
#'                tree_seeds = 1:25)
#'
#' orsf_vs(object)

orsf_vs <- function(object,
                    n_predictor_min = 3,
                    verbose_progress = FALSE){

 forest_seeds <- get_tree_seeds(object)

 if(is_empty(forest_seeds)){
  stop("tree_seeds not found in object. See details in ?orsf_vs",
       call. = FALSE)
 }

 forest_weights <- get_weights_user(object)

 if(is_empty(forest_weights))
  forest_weights <- NULL

 forest_outcomes <- get_names_y(object)
 formula <- stats::as.formula(
  paste( paste(forest_outcomes, collapse = ' + '), "~ ." )
 )

 forest_object <- object
 forest_data <- object$data
 forest_predictors <- get_names_x(object)
 n_predictors <- length(forest_predictors)

 oob_data <- data.table(
  n_predictors = seq(n_predictors),
  stat_value = rep(NA_real_, n_predictors),
  predictors_included = vector(mode = 'list', length = n_predictors),
  predictor_dropped = rep(NA_character_, n_predictors)
 )

 oob_last_stat_value <- get_last_oob_stat_value(forest_object)
 oob_worst_predictor <- get_last_vi(forest_object)

 oob_data[n_predictors,
          `:=`(n_predictors = n_predictors,
               stat_value = oob_last_stat_value,
               predictors_included = forest_predictors,
               predictor_dropped = oob_worst_predictor)]

 cols_kept <- c(
  forest_outcomes,
  setdiff(forest_predictors, oob_worst_predictor)
 )

 while(n_predictors > n_predictor_min){

  forest_data <- select_cols(forest_data, col_names = cols_kept)

  forest_object <- orsf(data = forest_data,
                        formula = formula,
                        control = get_control(object),
                        weights = forest_weights,
                        n_tree = get_n_tree(object),
                        n_split = get_n_split(object),
                        n_retry = get_n_retry(object),
                        leaf_min_events = get_leaf_min_events(object),
                        leaf_min_obs = get_leaf_min_obs(object),
                        split_min_events = get_split_min_events(object),
                        split_min_obs = get_split_min_obs(object),
                        split_min_stat = get_split_min_stat(object),
                        oobag_pred_type = get_oobag_pred_type(object),
                        oobag_pred_horizon = get_oobag_pred_horizon(object),
                        oobag_eval_every = get_n_tree(object),
                        oobag_fun = get_oobag_fun(object),
                        importance = get_importance(object),
                        group_factors = TRUE,
                        tree_seeds = forest_seeds,
                        attach_data = FALSE,
                        na_action = get_na_action(object),
                        verbose_progress = get_verbose_progress(object))

  forest_predictors <- get_names_x(forest_object)
  n_predictors <- length(forest_predictors)

  oob_last_stat_value <- get_last_oob_stat_value(forest_object)
  oob_worst_predictor <- get_last_vi(forest_object)

  oob_data[n_predictors,
           `:=`(n_predictors = n_predictors,
                stat_value = oob_last_stat_value,
                predictors_included = forest_predictors,
                predictor_dropped = oob_worst_predictor)]

  cols_kept <- c(
   forest_outcomes,
   setdiff(forest_predictors, get_last_vi(forest_object))
  )

 }

 collapse::na_omit(oob_data)

}


