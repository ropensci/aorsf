

#' Variable Interactions
#'
#' Use the variable interaction score described in Greenwell et al (2018).
#'   As this method can be computationally demanding, using `n_thread=0`
#'   can substantially reduce time needed to compute scores.
#'
#' @param object `r roxy_describe_ObliqueForest(TRUE)`
#'
#' @param predictors (*character*) a vector of length 2 or more with names
#'   of predictors used by `object`. All pairwise interactions between
#'   the predictors will be scored. If `NULL` (the default), all predictors
#'   are used.
#'
#' @param sep (*character*) how to separate the names of two predictors.
#'   The default value of `".."` returns names as `name1..name2`
#'
#' @inheritParams orsf
#'
#' @details
#'  The number of possible interactions grows exponentially based on the
#'  number of predictors. Some caution is warranted when using large predictor
#'  sets and it is recommended that you supply a specific vector of predictor
#'  names to assess rather than a global search. A good strategy is to use
#'  `n_tree = 5` to search all predictors, then pick the top 10 interactions,
#'  get the unique predictors from them, and re-run on just those predictors
#'  with more trees.
#'
#'
#' @return a data.table with variable interaction scores and
#'   partial dependence values.
#'
#' @export
#'
#' @references
#'
#' 1. `r cite("greenwell_2018")`
#'
#' @examples
#'
#' set.seed(329)
#'
#' data <- data.frame(
#'  x1 = rnorm(500),
#'  x2 = rnorm(500),
#'  x3 = rnorm(500)
#' )
#'
#' data$y = with(data, expr = x1 + x2 + x3 + 1/2*x1 * x2 + x2 * x3 + rnorm(500))
#'
#' forest <- orsf(data, y ~ ., n_tree = 5)
#'
#' orsf_vint(forest)

orsf_vint <- function(object,
                      predictors = NULL,
                      n_thread = NULL,
                      verbose_progress = NULL,
                      sep = '..'){

 check_arg_is(object, 'object', 'ObliqueForest')

 if(!is.null(predictors)){

  check_arg_type(arg_value = predictors,
                 arg_name = 'predictors',
                 expected_type = 'character')

 }

 pspec <- predictors %||% object$get_names_x()

 if(length(pspec) < 2)
  stop("at least two predictors are required.", call. = FALSE)

 class(pspec) <- c("pspec_intr", class(pspec))

 ptype <- switch(object$tree_type,
                 'survival' = 'mort',
                 'classification' = 'prob',
                 'regression' = 'mean')

 pd <-
  object$compute_dependence(
   pd_data = NULL,
   pred_spec = pspec,
   pred_horizon = NULL,
   pred_type = ptype,
   na_action = object$na_action,
   expand_grid = FALSE,
   prob_values = NULL,
   prob_labels = NULL,
   boundary_checks = FALSE,
   n_thread = n_thread %||% object$n_thread,
   verbose_progress = verbose_progress %||% object$verbose_progress,
   oobag = TRUE,
   type_output = "smry"
  )

 pd$id_intr <- paste(pd$var_1_name, pd$var_2_name, sep = sep)

 if(object$tree_type == 'classification'){
  pd[, mean := log(mean+0.01)]
  pd[, class := paste0(class, "._aorsf.split_")]
 }

 split_vars <- switch(object$tree_type,
                      "survival" = "id_intr",
                      "classification" = c("class", "id_intr"),
                      "regression" = "id_intr")

 pd_split <- split(pd, by = split_vars)

 # for cran
 . <- score <- var_1_value <- var_2_value <- pd_values <- NULL

 pd_scores <- vapply(
  pd_split,
  function(dt){
   collapse::fmean(
    c(
     collapse::fsd(dt[, .(vi = collapse::fsd(mean)), by = var_1_value][["vi"]]),
     collapse::fsd(dt[, .(vi = collapse::fsd(mean)), by = var_2_value][["vi"]])
    )
   )
  },
  double(1)
 )

 out <- data.table(interaction = names(pd_scores),
                   score = as.numeric(pd_scores))

 if(object$tree_type == 'classification'){

  out[, class := tstrsplit(interaction,
                           "\\.\\_aorsf\\.split\\_\\.",
                           keep = 1L)]

  out[, interaction := tstrsplit(interaction,
                                 "\\.\\_aorsf\\.split\\_\\.",
                                 keep = 2L)]

  out <- out[, .(score = mean(score)), by = c('interaction')]

 }

 out[, pd_values := pd_split]

 out[order(-score)]

}
