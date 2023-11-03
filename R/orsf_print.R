



#' Inspect your ORSF model
#'
#' Printing an ORSF model tells you:
#' - Linear combinations: How were these identified?
#' - N observations: Number of rows in training data
#' - N events: Number of events in training data
#' - N trees: Number of trees in the forest
#' - N predictors total: Total number of columns in the predictor matrix
#' - N predictors per node: Number of variables used in linear combinations
#' - Average leaves per tree: A proxy for the depth of your trees
#' - Min observations in leaf: See `leaf_min_obs` in [orsf]
#' - Min events in leaf: See `leaf_min_events` in [orsf]
#' - OOB stat value: Out-of-bag error after fitting all trees
#' - OOB stat type: How was out-of-bag error computed?
#' - Variable importance: How was variable importance computed?
#'
#' @param x (*orsf_fit*) an oblique random survival forest (ORSF; see [orsf]).
#'
#' @param ... `r roxy_dots()`
#'
#' @return `x`, invisibly.
#'
#' @export
#'
#' @examples
#'
#' object <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 5)
#'
#' print(object)
#'
print.orsf_fit <- function(x, ...){

 control            <- get_control(x)
 info_tree_type     <- get_tree_type(x)
 info_n_obs         <- get_n_obs(x)
 info_n_events      <- get_n_events(x) %||% "TBD"
 info_n_tree        <- get_n_tree(x)
 info_n_leaves_mean <- round_magnitude(get_n_leaves_mean(x)) %O% "TBD"
 info_max_time      <- get_max_time(x)
 info_mtry          <- get_mtry(x) %||% "TBD"
 info_names_x       <- get_names_x(x)
 info_leaf_min_obs  <- get_leaf_min_obs(x)
 info_vi            <- get_importance(x)
 info_leaf_min_events <- get_leaf_min_events(x)

 info_model_type <- switch(info_tree_type,
                           'survival' = "Cox regression",
                           'regression' = "Linear regression",
                           'classification' = "Logistic regression")

 info_lincomb_type <- switch(control$lincomb_type,
                             'glm'    = NULL,
                             'net'    = "Penalized ",
                             'custom' = "Custom user function")


 if(control$lincomb_type != 'custom'){

  if(control$lincomb_iter_max == 1){
   info_lincomb_type <- "Accelerated "
  }

  info_lincomb_type <- paste0(info_lincomb_type, info_model_type)

 }

 info_oobag_type <- info_oobag_stat <- 'none'

 if(contains_oobag(x)){

  info_oobag_type <- x$eval_oobag$stat_type

  info_oobag_stat <- "Not estimated"

  if(is_trained(x)){
   info_oobag_stat <- round_magnitude(last_value(x$eval_oobag$stat_values))
  }

 }


 forest_text <- paste("Oblique random", info_tree_type, "forest")

 header <- paste0('---------- ', forest_text, '\n')

 if(!is_trained(x)){
  header <- paste0('Untrained ', tolower(forest_text), '\n')
 }

 cat(header,
     paste0('     Linear combinations: ', info_lincomb_type    ),
     paste0('          N observations: ', info_n_obs           ),
     paste0('                N events: ', info_n_events        ),
     paste0('                 N trees: ', info_n_tree          ),
     paste0('      N predictors total: ', length(info_names_x) ),
     paste0('   N predictors per node: ', info_mtry            ),
     paste0(' Average leaves per tree: ', info_n_leaves_mean   ),
     paste0('Min observations in leaf: ', info_leaf_min_obs    ),
     paste0('      Min events in leaf: ', info_leaf_min_events ),
     paste0('          OOB stat value: ', info_oobag_stat      ),
     paste0('           OOB stat type: ', info_oobag_type      ),
     paste0('     Variable importance: ', info_vi              ),
     '\n-----------------------------------------',
     sep = '\n')

 invisible(x)

}
