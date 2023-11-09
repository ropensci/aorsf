



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
print.ObliqueForest <- function(x, ...){

 x$print()
 invisible(x)

}




