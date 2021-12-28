

#' Oblique Random Survival Forest (ORSF)
#'
#' The oblique random survival forest (RSF) is an extension of the RSF
#'   algorithm developed by Ishwaran et al and maintained in the
#'   `RandomForestSRC` package. The only difference between oblique
#'   RSF and Ishwaran's RSF is that oblique RSFs use linear combinations
#'   of input variables instead of using the input variable as-is when
#'   growing new nodes in survival decision trees. For more details on
#'   the oblique RSF, see Jaeger et al, 2019.
#'
#' @param data_train (_data.frame_) that will be used to grow the forest.
#'
#' @param formula (_formula_) a formula object, with the response on the left
#'   of a `~` operator, and the terms on the right. See details.
#'
#' @param n_tree (_integer_) the number of trees to grow
#'
#' @param n_split (_integer_) the number of cut-points assessed when splitting
#'  a node in decision trees.
#'
#' @param n_retry (_integer_) when a node can be split, but the current
#'  linear combination of inputs is unable to provide a valid split, `orsf`
#'  will try again with a new linear combination based on a different set
#'  of randomly selected predictors, up to `n_retry` times. When
#'  `n_retry = 0` (the default) the retry mechanic is not applied.
#'
#' @param mtry (_integer_) Number of variables randomly selected as candidates
#'   for splitting a node. The default is the smallest integer greater than
#'   the square root of the number of features. See details.
#'
#' @param leaf_min_events (_integer_) minimum number of events in a
#'   leaf node.
#'
#' @param leaf_min_obs (_integer_) minimum number of observations in a
#'   leaf node.
#'
#' @param oobag_pred (_logical_) if `TRUE` out-of-bag predictions are returned
#'   in the `aorsf` object.
#'
#' @param oobag_time (_numeric_) A numeric value indicating what time
#'   should be used for out-of-bag predictions.
#'
#' @param oobag_eval_every (_integer_) The out-of-bag performance of the
#'   ensemble will be checked every `oobag_eval_every` trees. So, if
#'   `oobag_eval_every = 10`, then out-of-bag performance is checked
#'   after growing the 10th tree, the 20th tree, and so on.
#'
#' @param importance (_logical_) if `TRUE`, variable importance will be
#'   computed using _negation_ importance. With negation importance,
#'   all coefficients for a given variable are multiplied by -1 and
#'   then the out-of-bag error for the forest is re-computed. The greater
#'   the degradation of the forest's error, the more important the variable.
#'
#' @param attach_data (_logical_) if `TRUE`, a copy of the training
#'   data will be attached to the output. This is helpful if you
#'   plan on using functions like [orsf_pd_summary] to interpret the fitted
#'   forest using its training data.
#'
#' @return an accelerated oblique RSF object (`aorsf`)
#'
#' @details
#'
#' This function is based on and highly similar to the `ORSF` function
#'   in the `obliqueRSF` R package. The primary difference is that this
#'   function runs about 200 times faster because it uses a simplified
#'   Newton Raphson scoring algorithm to identify linear combinations of
#'   inputs rather than performing penalized regression using routines in
#'   `glmnet`.The modified Newton Raphson scoring algorithm that this
#'   function applies is an adaptation of the C++ routine developed by
#'   Terry M. Therneau that fits Cox proportional hazards models
#'   (see [survival::coxph()]).
#'
#'
#' __What is an oblique decision tree?__
#'
#' Decision trees are developed by splitting a set of training data into two
#'  new subsets, with the goal of having more similarity within the new subsets
#'  than between them. This splitting process is repeated on the resulting
#'  subsets of data until a stopping criterion is met. When the new subsets of
#'  data are formed based on a single predictor, the decision tree is said to
#'  be axis-based because the splits of the data appear perpendicular to the
#'  axis of the predictor. When linear combinations of variables are used
#'  instead of a single variable, the tree is oblique because the splits of
#'  the data are neither parallel nor at a right angle to the axis
#'
#'
#' _Figure_ : Decision trees for classification with axis-based splitting
#'  (left) and oblique splitting (right). Cases are orange squares; controls
#'  are purple circles. Both trees partition the predictor space defined by
#'  variables X1 and X2, but the oblique splits do a better job of separating
#'  the two classes.
#'
#' \if{html}{\figure{tree_axis_v_oblique.png}{options: width=95\%}}
#'
#' __Some comments on inputs__
#'
#' _formula_: The response in `formula` can be a survival
#'   object as returned by the [survival::Surv] function,
#'   but can also just be the time and status variables.
#'   For example, `Surv(time, status) ~ .` works just like
#'   `time + status ~ .`. The only thing that can break this
#'   input is putting the variables in the wrong order, i.e.,
#'   writing `status + time ~ .` will make `orsf` assume your
#'   `status` variable is actually the `time` variable.
#'
#' _mtry_: The `mtry` parameter may be decreased while fitting an oblique
#'   RSF. Currently oblique RSF's are fitted by a Newton Raphson scoring
#'   algorithm that becomes unstable when the number of covariates is
#'   greater than or equal to the number of events. During the ORSF
#'   algorithm, mtry may be reduced temporarily to ensure there are at
#'   least 2 events per predictor variable.
#'
#' @references Jaeger BC, Long DL, Long DM, Sims M, Szychowski JM, Min YI,
#'   Mcclure LA, Howard G, Simon N. Oblique random survival forests.
#'   *The Annals of Applied Statistics*. 2019 Sep;13(3):1847-83.
#'   DOI: 10.1214/19-AOAS1261
#'
#' @export
#'
#' @examples
#'
#' fit <- orsf(pbc_orsf, formula = Surv(time, status) ~ . - id)
#'
#' print(fit)
#'
orsf <- function(data_train,
                 formula,
                 control = orsf_control_cph(),
                 n_tree = 500,
                 n_split = 5,
                 n_retry = 0,
                 mtry = NULL,
                 leaf_min_events = 1,
                 leaf_min_obs = 15,
                 oobag_pred = TRUE,
                 oobag_time = NULL,
                 oobag_eval_every = n_tree,
                 importance = FALSE,
                 attach_data = TRUE){


 orsf_type <- attr(control, 'type')

 switch(
  orsf_type,
  'cph' = {
   control_net <- orsf_control_net()
   control_cph <- control

  },
  'net' = {
   control_net <- control
   control_cph <- orsf_control_cph(do_scale = FALSE)
  })

 list2env(control_net, envir = environment())
 list2env(control_cph, envir = environment())

 # Run checks
 check_orsf_inputs(
  data_train = data_train,
  formula = formula,
  n_tree = n_tree,
  n_split = n_split,
  n_retry = n_retry,
  mtry = mtry,
  leaf_min_events = leaf_min_events,
  leaf_min_obs = leaf_min_obs,
  cph_method = cph_method,
  cph_eps = cph_eps,
  cph_iter_max = cph_iter_max,
  cph_pval_max = cph_pval_max,
  cph_do_scale = cph_do_scale,
  oobag_pred = oobag_pred,
  oobag_time = oobag_time,
  oobag_eval_every = oobag_eval_every,
  importance = importance,
  attach_data = attach_data
 )

 if(importance && !oobag_pred) oobag_pred <- TRUE # Should I add a warning?

 formula_terms <- suppressWarnings(stats::terms(formula, data=data_train))

 if(attr(formula_terms, 'response') == 0)
  stop("formula must have a response", call. = FALSE)

 if(length(attr(formula_terms, 'term.labels')) < 2)
  stop("formula must have at least 2 predictors", call. = FALSE)

 names_y_data <- all.vars(formula[[2]])

 if(length(names_y_data) != 2)
  stop("formula must have two variables (time & status) as the response",
       call. = FALSE)


 names_x_data <- attr(formula_terms, 'term.labels')

 names_not_found <- setdiff(c(names_y_data, names_x_data), names(data_train))

 if(!is_empty(names_not_found)){
  msg <- paste0(
   "variables in formula were not found in data_train: ",
   paste_collapse(names_not_found, last = ' and ')
  )
  stop(msg, call. = FALSE)
 }

 names_x_in_f <- rownames(attr(formula_terms, 'factors'))[-1]

 names_strange <- setdiff(names_x_in_f, names(data_train))

 if(!is_empty(names_strange)){
  msg <- paste0(
   "variables in formula were not found in data_train: ",
   paste_collapse(names_strange, last = ' and ')
  )
  warning(msg, call. = FALSE)
 }

 if(any(is.na(data_train[, c(names_y_data, names_x_data)]))){
  stop("Please remove missing values from data_train, or impute them.",
       call. = FALSE)
 }

 fctr_check(data_train, names_x_data)
 fctr_id_check(data_train, names_x_data)

 fi <- fctr_info(data_train, names_x_data)
 y  <- as.matrix(data_train[, names_y_data])
 x  <- as.matrix(one_hot(data_train, fi, names_x_data))

 types_x_data <- check_var_types(data_train,
                                 names_x_data,
                                 valid_types = c('numeric','integer',
                                                 'factor', 'ordered'))

 names_x_numeric <- grep(pattern = "^integer$|^numeric$",
                         x = types_x_data)

 numeric_bounds <- NULL

 if(!is_empty(names_x_numeric)){
  numeric_bounds <- sapply(data_train[, names_x_data[names_x_numeric] ],
                           FUN = stats::quantile,
                           probs = c(0.10, 0.25, 0.50, 0.75, 0.90))
 }


 if(is.null(mtry)){
  mtry <- ceiling(sqrt(ncol(x)))
 }

 if(is.null(net_df_target)) net_df_target <- mtry

 # TODO: expand
 if(net_df_target > mtry) stop("net_df_target > mtry")

 # Check the outcome variable


 check_arg_type(arg_value = y[, 2],
                arg_name = "status indicator",
                expected_type = 'numeric')

 check_arg_uni(arg_value = y[, 2],
               arg_name = "status indicator",
               expected_uni = c(0,1))

 n_events <- sum(y[, 2])

 check_arg_type(arg_value = y[, 1],
                arg_name = "time to event",
                expected_type = 'numeric')

 check_arg_gt(arg_value = y[, 1],
              arg_name = "time to event",
              bound = 0)

 # some additional checks that are dependent on the outcome variable

 check_arg_lteq(
  arg_value = mtry,
  arg_name = 'mtry',
  bound = ncol(x),
  append_to_msg = "(number of columns in the one-hot encoded x-matrix)"
 )

 check_arg_lteq(
  arg_value = leaf_min_events,
  arg_name = 'leaf_min_events',
  bound = round(n_events / 2),
  append_to_msg = "(number of events divided by 2)"
 )

 check_arg_lteq(
  arg_value = leaf_min_obs,
  arg_name = 'leaf_min_obs',
  bound = round(nrow(x) / 2),
  append_to_msg = "(number of observations divided by 2)"
 )




 if(!is.null(oobag_time)){

  if(oobag_time <= 0)

   stop("Out of bag prediction time (oobag_time) must be > 0",
        call. = FALSE)

 } else {

  # sneaky way to tell orsf.cpp to make its own oobag_time
  oobag_time <- 0

 }

 sorted <- order(y[, 1], -y[, 2])

 x_sort <- x[sorted, ]
 y_sort <- y[sorted, ]


 orsf_out <- orsf_fit(x                 = x_sort,
                      y                 = y_sort,
                      n_tree            = n_tree,
                      n_split_          = n_split,
                      mtry_             = mtry,
                      leaf_min_events_  = leaf_min_events,
                      leaf_min_obs_     = leaf_min_obs,
                      cph_method_       = switch(tolower(cph_method),
                                                 'breslow' = 0,
                                                 'efron'   = 1),
                      cph_eps_          = cph_eps,
                      cph_iter_max_     = cph_iter_max,
                      cph_pval_max_     = cph_pval_max,
                      cph_do_scale_     = cph_do_scale,
                      net_alpha_        = net_alpha,
                      net_df_target_    = net_df_target,
                      oobag_pred_       = oobag_pred,
                      oobag_time_       = oobag_time,
                      oobag_eval_every_ = oobag_eval_every,
                      oobag_importance_ = importance,
                      max_retry_        = n_retry,
                      penalized_cph     = penalized_cph,
                      type_             = switch(orsf_type,
                                                 'cph' = 'N',
                                                 'net' = 'P'))

 orsf_out$data_train <- if(attach_data) data_train else NULL

 if(importance){
  rownames(orsf_out$importance) <- colnames(x)
  orsf_out$importance <-
   rev(orsf_out$importance[order(orsf_out$importance), , drop=TRUE])
 }


 if(oobag_pred){

  # put the oob predictions into the same order as the training data.
  unsorted <- vector(mode = 'integer', length = length(sorted))
  for(i in seq_along(unsorted)) unsorted[ sorted[i] ] <- i

  orsf_out$surv_oobag <- orsf_out$surv_oobag[unsorted, , drop = FALSE]

 } else {

  # this would get added by orsf_fit if oobag_pred was TRUE
  orsf_out$time_pred <- stats::median(y[, 1])

 }

 n_leaves_mean <-
  mean(sapply(orsf_out$forest, function(t) nrow(t$leaf_node_index)))

 class(orsf_out) <- "aorsf"

 attr(orsf_out, 'mtry')            <- mtry
 attr(orsf_out, 'n_obs')           <- nrow(y_sort)
 attr(orsf_out, 'n_tree')          <- n_tree
 attr(orsf_out, 'names_y')         <- names_y_data
 attr(orsf_out, "names_x")         <- names_x_data
 attr(orsf_out, "names_x_onehot")  <- colnames(x)
 attr(orsf_out, "types_x")         <- types_x_data
 attr(orsf_out, 'n_events')        <- n_events
 attr(orsf_out, 'max_time')        <- y_sort[nrow(y_sort), 1]
 attr(orsf_out, "fctr_info")       <- fi
 attr(orsf_out, 'n_leaves_mean')   <- n_leaves_mean
 attr(orsf_out, 'n_split')         <- n_split
 attr(orsf_out, 'leaf_min_events') <- leaf_min_events
 attr(orsf_out, 'leaf_min_obs')    <- leaf_min_obs
 attr(orsf_out, 'numeric_bounds')  <- numeric_bounds

 orsf_out


}

check_orsf_inputs <- function(data_train,
                              formula,
                              n_tree,
                              n_split,
                              n_retry,
                              mtry,
                              leaf_min_events,
                              leaf_min_obs,
                              cph_method,
                              cph_eps,
                              cph_iter_max,
                              cph_pval_max,
                              cph_do_scale,
                              oobag_pred,
                              oobag_time,
                              oobag_eval_every,
                              importance,
                              attach_data){

 if(!is.null(data_train)){

  check_arg_is(arg_value = data_train,
               arg_name = 'data_train',
               expected_class = 'data.frame')

 }

 if(!is.null(formula)){

  check_arg_is(arg_value = formula,
               arg_name = 'formula',
               expected_class = 'formula')

  if(length(formula) != 3){
   stop("formula must be two sided, i.e. left side ~ right side",
        call. = FALSE)
  }

  formula_deparsed <- deparse(formula[[3]])

  for( symbol in c("*", "^", ":", "(", ")", "["," ]", "|", "%") ){

   if(grepl(symbol, formula_deparsed, fixed = TRUE)){

    stop("unrecognized symbol in formula: ", symbol,
         "\norsf recognizes '+', '-', and '.' symbols.",
         call. = FALSE)

   }

  }

 }

 if(!is.null(n_tree)){

  check_arg_type(arg_value = n_tree,
                 arg_name = 'n_tree',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = n_tree,
                       arg_name = 'n_tree')

  check_arg_gteq(arg_value = n_tree,
                 arg_name = 'n_tree',
                 bound = 1)

  check_arg_length(arg_value = n_tree,
                   arg_name = 'n_tree',
                   expected_length = 1)

 }

 if(!is.null(n_split)){

  check_arg_type(arg_value = n_split,
                 arg_name = 'n_split',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = n_split,
                       arg_name = 'n_split')

  check_arg_gteq(arg_value = n_split,
                 arg_name = 'n_split',
                 bound = 1)

  check_arg_length(arg_value = n_split,
                   arg_name = 'n_split',
                   expected_length = 1)

 }

 if(!is.null(n_retry)){

  check_arg_type(arg_value = n_retry,
                 arg_name = 'n_retry',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = n_retry,
                       arg_name = 'n_retry')

  check_arg_gteq(arg_value = n_retry,
                 arg_name = 'n_retry',
                 bound = 0)

  check_arg_length(arg_value = n_retry,
                   arg_name = 'n_retry',
                   expected_length = 1)

 }

 if(!is.null(mtry)){

  check_arg_type(arg_value = mtry,
                 arg_name = 'mtry',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_name = 'mtry',
                       arg_value = mtry)

  check_arg_gteq(arg_name = 'mtry',
                 arg_value = mtry,
                 bound = 2)

  check_arg_length(arg_name = 'mtry',
                   arg_value = mtry,
                   expected_length = 1)

 }

 if(!is.null(leaf_min_events)){

  check_arg_type(arg_value = leaf_min_events,
                 arg_name = 'leaf_min_events',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = leaf_min_events,
                       arg_name = 'leaf_min_events')

  check_arg_gteq(arg_value = leaf_min_events,
                 arg_name = 'leaf_min_events',
                 bound = 1)

  check_arg_length(arg_value = leaf_min_events,
                   arg_name = 'leaf_min_events',
                   expected_length = 1)
 }

 if(!is.null(leaf_min_obs)){

  check_arg_type(arg_value = leaf_min_obs,
                 arg_name = 'leaf_min_obs',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = leaf_min_obs,
                       arg_name = 'leaf_min_obs')

  check_arg_gteq(arg_value = leaf_min_obs,
                 arg_name = 'leaf_min_obs',
                 bound = 1)

  check_arg_length(arg_value = leaf_min_obs,
                   arg_name = 'leaf_min_obs',
                   expected_length = 1)

 }

 if(!is.null(cph_method)){

  check_arg_type(arg_value = cph_method,
                 arg_name = 'cph_method',
                 expected_type = 'character')

  check_arg_is_valid(arg_value = cph_method,
                     arg_name = 'cph_method',
                     valid_options = c("breslow", "efron"))
 }

 if(!is.null(cph_eps)){

  check_arg_type(arg_value = cph_eps,
                 arg_name = 'cph_eps',
                 expected_type = 'numeric')

  check_arg_gt(arg_value = cph_eps,
               arg_name = 'cph_eps',
               bound = 0)

  check_arg_length(arg_value = cph_eps,
                   arg_name = 'cph_eps',
                   expected_length = 1)

 }

 if(!is.null(cph_iter_max)){

  check_arg_type(arg_value = cph_iter_max,
                 arg_name = 'cph_iter_max',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = cph_iter_max,
                       arg_name = 'cph_iter_max')

  check_arg_gteq(arg_value = cph_iter_max,
                 arg_name = 'cph_iter_max',
                 bound = 1)

  check_arg_length(arg_value = cph_iter_max,
                   arg_name = 'cph_iter_max',
                   expected_length = 1)

 }

 if(!is.null(cph_pval_max)){

  check_arg_type(arg_value = cph_pval_max,
                 arg_name = 'cph_pval_max',
                 expected_type = 'numeric')

  check_arg_gt(arg_value = cph_pval_max,
               arg_name = 'cph_pval_max',
               bound = 0)

  check_arg_lteq(arg_value = cph_pval_max,
                 arg_name = 'cph_pval_max',
                 bound = 1)

  check_arg_length(arg_value = cph_pval_max,
                   arg_name = 'cph_pval_max',
                   expected_length = 1)

 }

 if(!is.null(cph_do_scale)){

  check_arg_type(arg_value = cph_do_scale,
                 arg_name = 'cph_do_scale',
                 expected_type = 'logical')

  check_arg_length(arg_value = cph_do_scale,
                   arg_name = 'cph_do_scale',
                   expected_length = 1)

 }

 if(!is.null(oobag_pred)){

  check_arg_type(arg_value = oobag_pred,
                 arg_name = 'oobag_pred',
                 expected_type = 'logical')

  check_arg_length(arg_value = oobag_pred,
                   arg_name = 'oobag_pred',
                   expected_length = 1)

 }

 if(!is.null(oobag_time)){

  check_arg_type(arg_value = oobag_time,
                 arg_name = 'oobag_time',
                 expected_type = 'numeric')

  check_arg_length(arg_value = oobag_time,
                   arg_name = 'oobag_time',
                   expected_length = 1)

  check_arg_gt(arg_value = oobag_time,
               arg_name = 'oobag_time',
               bound = 0)

 }


 if(!is.null(oobag_eval_every)){

  check_arg_type(arg_value = oobag_eval_every,
                 arg_name = 'oobag_eval_every',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = oobag_eval_every,
                       arg_name = 'oobag_eval_every')

  check_arg_gteq(arg_value = oobag_eval_every,
                 arg_name = 'oobag_eval_every',
                 bound = 1)

  check_arg_lteq(arg_value = oobag_eval_every,
                 arg_name = 'oobag_eval_every',
                 bound = n_tree)

  check_arg_length(arg_value = oobag_eval_every,
                   arg_name = 'oobag_eval_every',
                   expected_length = 1)

 }

 if(!is.null(attach_data)){

  check_arg_type(arg_value = attach_data,
                 arg_name = 'attach_data',
                 expected_type = 'logical')

  check_arg_length(arg_value = attach_data,
                   arg_name = 'attach_data',
                   expected_length = 1)

 }


 if(!cph_do_scale && cph_iter_max > 1){
  stop("cph_do_scale must be TRUE when cph_iter_max > 1",
       call. = FALSE)
 }


}









