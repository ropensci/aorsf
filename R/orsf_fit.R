

#' Oblique Random Survival Forest (ORSF)
#'
#' @srrstats {G1.4} *documented with Roxygen*
#'
#' @srrstats {G1.1} *aorsf is an improvement of the ORSF algorithm implemented in obliqueRSF, which was an extension of Hemant Ishwaran's random survival forest.*
#'
#' @srrstats {G1.3} *linear combinations of inputs defined.*
#'
#' The oblique random survival forest (ORSF) is an extension of the RSF
#'   algorithm developed by Ishwaran et al and maintained in the
#'   `RandomForestSRC` package. The difference between ORSF and RSF is
#'   that ORSF uses linear combinations of input variables whereas RSF
#'   uses a single variable when growing new nodes in survival decision trees.
#'   A linear combination is an expression constructed from a set of terms
#'   by multiplying each term by a constant and adding the results (e.g.
#'   a linear combination of x and y would be any expression of the form
#'   ax + by, where a and b are constants). For more details on the ORSF
#'   algorithm, see Jaeger et al, 2019. The `orsf()` function implements a
#'   novel algorithm that speeds up the ORSF algorithm described by Jaeger
#'   et al (see details).
#'
#' @param data_train (_data.frame_) that will be used to grow the forest.
#'
#' @param formula (_formula_) a formula object, with the response on the left
#'   of a `~` operator, and the terms on the right. See details.
#'
#' @param control An `aorsf_control` object, created with [orsf_control_net]
#'  or [orsf_control_cph]. Default is `orsf_control_cph()`.
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
#' @param split_min_events (_integer_) minimum number of events required
#'   to split a node.
#'
#' @param split_min_obs (_integer_) minimum number of observations required
#'   to split a node.
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
#'   function runs about 500 times faster because it uses a simplified
#'   Newton Raphson scoring algorithm to identify linear combinations of
#'   inputs rather than performing penalized regression using routines in
#'   `glmnet`.The modified Newton Raphson scoring algorithm that this
#'   function applies is an adaptation of the C++ routine developed by
#'   Terry M. Therneau that fits Cox proportional hazards models
#'   (see [survival::coxph()]).
#'
#'
#' @srrstats {G1.3} *define oblique and axis based decision trees*
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
#' @srrstats {G1.3} *clarify the term 'random forest'*
#'
#' __What is a random forest?__
#'
#' Random forests are collections of de-correlated decision trees. Predictions from each tree are aggregated to make an ensemble prediction for the forest. For more details, see Breiman at el, 2001.
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
#' _mtry_: The `mtry` parameter may be temporarily reduced to ensure there
#'   are at least 2 events per predictor variable. This occurs when using
#'   [orsf_control_cph] because coefficients in the Newton Raphson scoring
#'   algorithm may become unstable when the number of covariates is
#'   greater than or equal to the number of events. This reduction does not
#'   occur when using [orsf_control_net].
#'
#'
#' @srrstats {G1.0} *Jaeger et al describes the ORSF algorithm that aorsf is based on. Note: aorsf uses a different approach to create linear combinations of inputs for speed reasons, but orsf_control_net() allows users to make ensembles that are very similar to obliqueRSF::ORSF().*
#'
#' @references
#'
#' Breiman L. Random forests. *Machine learning*. 2001 Oct;45(1):5-32.
#'   DOI: 10.1023/A:1010933404324
#'
#' Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS. Random survival forests.
#'   *Annals of applied statistics*. 2008 Sep;2(3):841-60.
#'   DOI: 10.1214/08-AOAS169
#'
#' Jaeger BC, Long DL, Long DM, Sims M, Szychowski JM, Min YI,
#'   Mcclure LA, Howard G, Simon N. Oblique random survival forests.
#'   *Annals of applied statistics*. 2019 Sep;13(3):1847-83.
#'   DOI: 10.1214/19-AOAS1261
#'
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
                 leaf_min_obs = 5,
                 split_min_events = 5,
                 split_min_obs = 10,
                 oobag_pred = TRUE,
                 oobag_time = NULL,
                 oobag_eval_every = n_tree,
                 importance = FALSE,
                 attach_data = TRUE){


 # Run checks
 check_orsf_inputs(
  data_train = data_train,
  formula = formula,
  control = control,
  n_tree = n_tree,
  n_split = n_split,
  n_retry = n_retry,
  mtry = mtry,
  leaf_min_events = leaf_min_events,
  leaf_min_obs = leaf_min_obs,
  split_min_events = split_min_events,
  split_min_obs = split_min_obs,
  oobag_pred = oobag_pred,
  oobag_time = oobag_time,
  oobag_eval_every = oobag_eval_every,
  importance = importance,
  attach_data = attach_data
 )

 orsf_type <- attr(control, 'type')

 switch(
  orsf_type,

  'cph' = {
   control_net <- orsf_control_net()
   control_cph <- control
   f_beta      <- function(x) x

  },

  'net' = {

   if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop(
     "Package \"glmnet\" must be installed to use",
     " orsf_control_net() with orsf().",
     call. = FALSE
    )
   }

   control_net <- control
   control_cph <- orsf_control_cph(do_scale = FALSE)
   f_beta      <- penalized_cph
  }

 )

 cph_method = control_cph$cph_method
 cph_eps = control_cph$cph_eps
 cph_iter_max = control_cph$cph_iter_max
 cph_pval_max = control_cph$cph_pval_max
 cph_do_scale = control_cph$cph_do_scale

 net_alpha = control_net$net_alpha
 net_df_target = control_net$net_df_target

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
 y  <- as.matrix(select_cols(data_train, names_y_data))
 x  <- as.matrix(ref_code(data_train, fi, names_x_data))

 types_x_data <- check_var_types(data_train,
                                 names_x_data,
                                 valid_types = c('numeric','integer',
                                                 'factor', 'ordered'))

 names_x_numeric <- grep(pattern = "^integer$|^numeric$",
                         x = types_x_data)

 numeric_bounds <- NULL

 if(!is_empty(names_x_numeric)){
  numeric_bounds <-
   sapply(select_cols(data_train, names_x_data[names_x_numeric]),
          FUN = stats::quantile,
          probs = c(0.10, 0.25, 0.50, 0.75, 0.90))
 }


 if(is.null(mtry)){
  mtry <- ceiling(sqrt(ncol(x)))
 }

 if(is.null(net_df_target)) net_df_target <- mtry

 # TODO: warn instead?
 if(net_df_target > mtry)
  stop("net_df_target = ", net_df_target,
       " must be <= mtry, which is ", mtry,
       call. = FALSE)

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

 check_arg_lt(arg_value = split_min_events,
              arg_name = "split_min_events",
              bound = n_events,
              append_to_msg = "(number of events)")

 check_arg_lt(arg_value = split_min_obs,
              arg_name = "split_min_obs",
              bound = nrow(x),
              append_to_msg = "(number of observations)")

 if(!is.null(oobag_time)){

  if(oobag_time <= 0)

   stop("Out of bag prediction time (oobag_time) must be > 0",
        call. = FALSE)

 } else {

  # sneaky way to tell orsf.cpp to make its own oobag_time
  oobag_time <- 0

 }


 sorted <- order(y[, 1],  # order this way for risk sets
                 -y[, 2]) # order this way for oob C-statistic.

 x_sort <- x[sorted, ]
 y_sort <- y[sorted, ]

 orsf_out <- orsf_fit(x                 = x_sort,
                      y                 = y_sort,
                      n_tree            = n_tree,
                      n_split_          = n_split,
                      mtry_             = mtry,
                      leaf_min_events_  = leaf_min_events,
                      leaf_min_obs_     = leaf_min_obs,
                      split_min_events_ = split_min_events,
                      split_min_obs_    = split_min_obs,
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
                      f_beta            = f_beta,
                      type_             = switch(orsf_type,
                                                 'cph' = 'C',
                                                 'net' = 'N'))

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
 attr(orsf_out, "names_x_ref")  <- colnames(x)
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












