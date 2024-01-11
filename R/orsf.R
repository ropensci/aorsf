
#' Oblique Random Forests
#'
#' Grow or specify an oblique random forest. While the name `orsf()`
#'  implies that this function only works for survival forests,
#'  it can be used for classification, regression, or survival
#'  forests.
#'
#' @param data a `r roxy_data_allowed()` that contains the
#'  relevant variables.
#'
#' @param formula (*formula*) Two sided formula with a single outcome.
#'   The terms on the right are names of predictor variables, and the
#'   symbol '.' may be used to indicate all variables in the data
#'   except the response. The symbol '-' may also be used to indicate
#'   removal of a predictor. Details on the response vary depending
#'   on forest type:
#'
#'   - *Classification*: The response should be a single variable,
#'   and that variable should have type `factor` in `data`.
#'
#'   - *Regression*: The response should be a single variable, and
#'   that variable should have typee `double` or `integer` with at
#'   least 10 unique numeric values in `data`.
#'
#'   - *Survival*: The response should include a time variable,
#'   followed by a status variable, and may be written inside a
#'   call to [Surv][survival::Surv] (see examples).
#'
#'
#' @param control (*orsf_control*) An object returned from one of the
#'  `orsf_control` functions: [orsf_control_survival],
#'  [orsf_control_classification], and [orsf_control_regression]. If
#'  `NULL` (the default) will use an accelerated control, which is the
#'   fastest available option. For survival and classification, this is
#'   Cox and Logistic regression with 1 iteration, and for regression
#'   it is ordinary least squares.
#'
#' @param weights (*numeric vector*) Optional. If given, this input should
#'   have length equal to `nrow(data)` for complete or imputed data and should
#'   have length equal to `nrow(na.omit(data))` if `na_action` is `"omit"`.
#'   As the weights vector is used to count observations and events prior to
#'   growing a node for a tree, `orsf()` scales `weights` so that
#'   `sum(weights) == nrow(data)`. This helps to make tree depth consistent
#'   between weighted and un-weighted fits.
#'
#' @param n_tree (*integer*) the number of trees to grow.
#' Default is `n_tree = 500.`
#'
#' @param n_split (*integer*) the number of cut-points assessed when splitting
#'  a node in decision trees. Default is `n_split = 5`.
#'
#' @param n_retry (*integer*) when a node is splittable, but the current
#'  linear combination of inputs is unable to provide a valid split, `orsf`
#'  will try again with a new linear combination based on a different set
#'  of randomly selected predictors, up to `n_retry` times. Default is
#'  `n_retry = 3`. Set `n_retry = 0` to prevent any retries.
#'
#' @param n_thread `r roxy_n_thread_header("growing trees, computing predictions, and computing importance")`
#'
#' @param mtry (*integer*) Number of predictors randomly included as candidates
#'   for splitting a node. The default is the smallest integer greater than
#'   the square root of the number of total predictors, i.e.,
#'   `mtry = ceiling(sqrt(number of predictors))`
#'
#' @param sample_with_replacement (*logical*) If `TRUE` (the default),
#'   observations are sampled with replacement when an in-bag sample
#'   is created for a decision tree. If `FALSE`, observations are
#'   sampled without replacement and each tree will have an in-bag sample
#'   containing `sample_fraction`% of the original sample.
#'
#' @param sample_fraction (*double*) the proportion of observations that
#'   each trees' in-bag sample will contain, relative to the number of
#'   rows in `data`. Only used if `sample_with_replacement` is `FALSE`.
#'   Default value is 0.632.
#'
#' @param leaf_min_events (*integer*) This input is only relevant for
#'   survival analysis, and specifies the minimum number of events in a
#'   leaf node. Default is `leaf_min_events = 1`
#'
#' @param leaf_min_obs (*integer*) minimum number of observations in a
#'   leaf node. Default is `leaf_min_obs = 5`.
#'
#' @param split_rule (*character*) how to assess the quality of a potential
#'   splitting rule for a node. Valid options for survival are:
#'
#'   - 'logrank' : a log-rank test statistic (default).
#'   - 'cstat'   : Harrell's concordance statistic.
#'
#'   For classification, valid options are:
#'
#'   - 'gini'  : gini impurity (default)
#'   - 'cstat' : area underneath the ROC curve (AUC-ROC)
#'
#'   For regression, valid options are:
#'
#'   - 'variance' : variance reduction (default)
#'
#' @param split_min_events (*integer*) minimum number of events required
#'   in a node to consider splitting it. Default is `split_min_events = 5`.
#'   This input is only relevant for survival trees.
#'
#' @param split_min_obs (*integer*) minimum number of observations required
#'   in a node to consider splitting it. Default is `split_min_obs = 10`.
#'
#' @param split_min_stat (double) minimum test statistic required to split
#'   a node. If no splits are found with a statistic exceeding `split_min_stat`,
#'   the given node either becomes a leaf or a retry occurs (up to `n_retry`
#'   retries). Defaults are
#'
#'   - 3.84 if `split_rule = 'logrank'`
#'   - 0.55 if `split_rule = 'cstat'` (see first note below)
#'   - 0.00 if `split_rule = 'gini'` (see second note below)
#'   - 0.00 if `split_rule = 'variance'`
#'
#'   **Note 1** For C-statistic splitting, if C is < 0.50, we consider the statistic
#'   value to be 1 - C to allow for good 'anti-predictive' splits. So,
#'   if a C-statistic is initially computed as 0.1, it will be considered
#'   as 1 - 0.10 = 0.90.
#'
#'   **Note 2** For Gini impurity, a value of 0 and 1 usually indicate the best and
#'   worst possible scores, respectively. To make things simple and to avoid
#'   introducing a `split_max_stat` input, we flip the values of Gini
#'   impurity so that 1 and 0 indicate the best and worst possible scores,
#'   respectively.
#'
#' @param oobag_pred_type (*character*) The type of out-of-bag predictions
#'   to compute while fitting the ensemble. Valid options for any tree type:
#'
#'   - 'none' : don't compute out-of-bag predictions
#'   - 'leaf' : the ID of the predicted leaf is returned for each tree
#'
#'   Valid options for survival:
#'
#'   - 'risk' : probability of event occurring at or before
#'              `oobag_pred_horizon` (default).
#'   - 'surv' : 1 - risk.
#'   - 'chf'  : cumulative hazard function at `oobag_pred_horizon`.
#'   - 'mort' : mortality, i.e., the number of events expected if all
#'              observations in the training data were identical to a
#'              given observation.
#'
#'   Valid options for classification:
#'
#'   - 'prob'  : probability of each class (default)
#'   - 'class' : class (i.e., which.max(prob))
#'
#'   Valid options for regression:
#'
#'   - 'mean' : mean value (default)
#'
#' @param oobag_pred_horizon (*numeric*) A numeric value indicating what time
#'   should be used for out-of-bag predictions. Default is the median
#'   of the observed times, i.e., `oobag_pred_horizon = median(time)`.
#'   This input is only relevant for survival trees that have prediction
#'   type of 'risk', 'surv', or 'chf'.
#'
#' @param oobag_eval_every (*integer*) The out-of-bag performance of the
#'   ensemble will be checked every `oobag_eval_every` trees. So, if
#'   `oobag_eval_every = 10`, then out-of-bag performance is checked
#'   after growing the 10th tree, the 20th tree, and so on. Default
#'   is `oobag_eval_every = n_tree`.
#'
#' @param oobag_fun `r roxy_oobag_fun_header()` every `oobag_eval_every`
#'  trees. `r roxy_oobag_fun_default()`
#'
#'  `r roxy_oobag_fun_user()`
#'   - `r roxy_oobag_fun_inputs()`
#'   - `r roxy_oobag_fun_ymat()`
#'   - `r roxy_oobag_fun_svec()`
#'   - `r roxy_oobag_fun_return()`
#'
#' For more details, see the out-of-bag [vignette](https://docs.ropensci.org/aorsf/articles/oobag.html#user-supplied-out-of-bag-evaluation-functions).
#'
#' @param importance `r roxy_importance_header()`
#' - `r roxy_importance_none()`
#' - `r roxy_importance_anova()`
#' - `r roxy_importance_negate()`
#' - `r roxy_importance_permute()`
#'
#' For details on these methods, see [orsf_vi].
#'
#' @param importance_max_pvalue (*double*) Only relevant if `importance`
#'   is `"anova"`. The maximum p-value that will register as a positive
#'   case when counting the number of times a variable was found to be
#'   'significant' during tree growth. Default is 0.01, as recommended
#'   by Menze et al.
#'
#' @param group_factors (*logical*) Only relevant if variable importance is
#'   being estimated. `r roxy_group_factors()`
#'
#' @param tree_seeds (*integer vector*) Optional. if specified, random seeds
#'   will be set using the values in `tree_seeds[i]`  before growing tree `i`.
#'   Two forests grown with the same number of trees and the same seeds will
#'   have the exact same out-of-bag samples, making out-of-bag error
#'   estimates of the forests more comparable. If `NULL` (the default),
#'   seeds are picked at random.
#'
#' @param attach_data (*logical*) if `TRUE`, a copy of the training
#'   data will be attached to the output. This is required if you
#'   plan on using functions like [orsf_pd_oob] or [orsf_summarize_uni]
#'   to interpret the forest using its training data. Default is `TRUE`.
#'
#' @param no_fit (*logical*) if `TRUE`, model fitting steps are defined and
#'   saved, but training is not initiated. The object returned can be
#'   directly submitted to `orsf_train()` so long as `attach_data` is `TRUE`.
#'
#' @param na_action `r roxy_na_action_header("data")`
#'
#'   - `r roxy_na_action_fail("data")`
#'   - `r roxy_na_action_omit("data")`
#'   - `r roxy_na_action_impute_meanmode("data")`.
#'
#' @param verbose_progress (*logical*) if `TRUE`, progress messages are
#'   printed in the console. If `FALSE` (the default), nothing is printed.
#'
#' @param ... `r roxy_dots()`
#'
#' @param object an untrained 'aorsf' object, created by setting
#'   `no_fit = TRUE` in `orsf()`.
#'
#' @return an *obliqueForest* object
#'
#' @details
#'
#' Why isn't this function called `orf()`? In its earlier versions, the
#' `aorsf` package was exclusively for *o*blique *r*andom *s*urvival *f*orests.
#'
#' **formula for survival oblique RFs**:
#'
#' - The response in `formula` can be a survival
#'   object as returned by the [Surv][survival::Surv] function,
#'   but can also just be the time and status variables. I.e.,
#'   `Surv(time, status) ~ .` works and `time + status ~ .` works
#'
#' - The response can also be a survival object stored in `data`.
#'   For example, `y ~ .` is a valid formula if `data$y` inherits
#'   from the `Surv` class.
#'
#' **mtry**:
#'
#' The `mtry` parameter may be temporarily reduced to ensure that linear
#'   models used to find combinations of predictors remain stable. This occurs
#'   because coefficients in linear model fitting algorithms may become infinite
#'   if the number of predictors exceeds the number of observations.
#'
#' **oobag_fun**:
#'
#' If `oobag_fun` is specified, it will be used in to compute negation
#'  importance or permutation importance, but it will not have any role
#'  for ANOVA importance.
#'
#' **n_thread**:
#'
#' If an R function is to be called from C++ (i.e., user-supplied function to
#'  compute out-of-bag error or identify linear combinations of variables),
#'  `n_thread` will automatically be set to 1 because attempting to run R
#'  functions in multiple threads will cause the R session to crash.
#'
#' @section What is an oblique decision tree?:
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
#' _Figure_ : Decision trees for classification with axis-based splitting
#'  (left) and oblique splitting (right). Cases are orange squares; controls
#'  are purple circles. Both trees partition the predictor space defined by
#'  variables X1 and X2, but the oblique splits do a better job of separating
#'  the two classes.
#'
#' \if{html}{\figure{tree_axis_v_oblique.png}{options: width=95\%}}
#'
#'
#' @section What is a random forest?:
#'
#' Random forests are collections of de-correlated decision trees.
#'   Predictions from each tree are aggregated to make an ensemble
#'   prediction for the forest. For more details, see Breiman at el, 2001.
#'
#' @section Training, out-of-bag error, and testing:
#'
#' In random forests, each tree is grown with a bootstrapped version of
#'   the training set. Because bootstrap samples are selected with replacement,
#'   each bootstrapped training set contains about two-thirds of instances in
#'   the original training set. The 'out-of-bag' data are instances that are
#'   _not_ in the bootstrapped training set. Each tree in the random forest
#'   can make predictions for its out-of-bag data, and the out-of-bag
#'   predictions can be aggregated to make an ensemble out-of-bag prediction.
#'   Since the out-of-bag data are not used to grow the tree, the accuracy of
#'   the ensemble out-of-bag predictions approximate the generalization error
#'   of the random forest. Generalization error refers to the error of a
#'   random forest's predictions when it is applied to predict outcomes for
#'   data that were not used to train it, i.e., testing data.
#'
#' @references
#'
#' 1. `r cite("harrell_1982")`
#' 1. `r cite("breiman_2001")`
#' 1. `r cite("ishwaran_2008")`
#' 1. `r cite("menze_2011")`
#' 1. `r cite("jaeger_2019")`
#' 1. `r cite("jaeger_2022")`
#'
#' @includeRmd Rmd/orsf_examples.Rmd
#'
#' @export
#'



orsf <- function(data,
                 formula,
                 control = NULL,
                 weights = NULL,
                 n_tree = 500,
                 n_split = 5,
                 n_retry = 3,
                 n_thread = 0,
                 mtry = NULL,
                 sample_with_replacement = TRUE,
                 sample_fraction = 0.632,
                 leaf_min_events = 1,
                 leaf_min_obs = 5,
                 split_rule = NULL,
                 split_min_events = 5,
                 split_min_obs = 10,
                 split_min_stat = NULL,
                 oobag_pred_type = NULL,
                 oobag_pred_horizon = NULL,
                 oobag_eval_every = NULL,
                 oobag_fun = NULL,
                 importance = 'anova',
                 importance_max_pvalue = 0.01,
                 group_factors = TRUE,
                 tree_seeds = NULL,
                 attach_data = TRUE,
                 no_fit = FALSE,
                 na_action = 'fail',
                 verbose_progress = FALSE,
                 ...){

 check_dots(list(...), .f = orsf)

 check_arg_type(arg_value = attach_data,
                arg_name = 'attach_data',
                expected_type = 'logical')

 check_arg_length(arg_value = attach_data,
                  arg_name = 'attach_data',
                  expected_length = 1)

 check_arg_type(arg_value = no_fit,
                arg_name = 'no_fit',
                expected_type = 'logical')

 check_arg_length(arg_value = no_fit,
                  arg_name = 'no_fit',
                  expected_length = 1)

 args <- list(data = orsf_data_prep(data),
              formula = formula,
              control = control,
              weights = weights,
              n_tree = n_tree,
              n_split = n_split,
              n_retry = n_retry,
              n_thread = n_thread,
              mtry = mtry,
              sample_with_replacement = sample_with_replacement,
              sample_fraction = sample_fraction,
              leaf_min_events = leaf_min_events,
              leaf_min_obs = leaf_min_obs,
              split_rule = split_rule,
              split_min_events = split_min_events,
              split_min_obs = split_min_obs,
              split_min_stat = split_min_stat,
              pred_type = oobag_pred_type,
              oobag_pred_horizon = oobag_pred_horizon,
              oobag_eval_every = oobag_eval_every,
              oobag_fun = oobag_fun,
              importance_type = importance,
              importance_max_pvalue = importance_max_pvalue,
              importance_group_factors = group_factors,
              tree_seeds = tree_seeds,
              na_action = na_action,
              verbose_progress = verbose_progress)

 if(length(formula) < 3) stop("formula must be two sided", call. = FALSE)

 if(is.null(control)){
  tree_type <- infer_tree_type(all.vars(formula[[2]]), args$data)
 } else {
  tree_type <- control$tree_type
 }


 object <- switch(
  tree_type,
  'survival' = do.call(ObliqueForestSurvival$new, args = args),
  'classification' = do.call(ObliqueForestClassification$new, args = args),
  'regression' = do.call(ObliqueForestRegression$new, args = args)
 )

 if(no_fit) return(object)

 orsf_train(object, attach_data = attach_data)

}

#' @rdname orsf
#' @export
orsf_train <- function(object, attach_data = TRUE){

 object$train()

 if(!attach_data) object$data <- NULL

 invisible(object)

}

#' Estimate training time
#'
#' @param object an untrained `aorsf` object
#'
#' @param n_tree_subset (*integer*)  how many trees should be fit in order
#'   to estimate the time needed to train `object`. The default value is 50,
#'   as this usually gives a good enough approximation.
#'
#' @return a [difftime] object.
#'
#' @export
#'
#' @examples
#'
#' # specify but do not train the model by setting no_fit = TRUE.
#' object <- orsf(pbc_orsf, Surv(time, status) ~ . - id,
#'                n_tree = 10, no_fit = TRUE)
#'
#' # approximate the time it will take to grow 500 trees
#' time_estimated <- orsf_time_to_train(object)
#'
#' print(time_estimated)
#'
#' # let's see how close the approximation was
#' time_true_start <- Sys.time()
#' orsf_train(object)
#' time_true_stop <- Sys.time()
#'
#' time_true <- time_true_stop - time_true_start
#'
#' print(time_true)
#'
#' # error
#' abs(time_true - time_estimated)
#'

orsf_time_to_train <- function(object, n_tree_subset = 50){

 n_tree_original <- object$n_tree

 time_train_start <- Sys.time()

 object$train(n_tree = n_tree_subset)

 time_train_stop <- Sys.time()

 object$untrain()

 object$set_field(n_tree = n_tree_original)

 difftime(time_train_stop, time_train_start) * object$n_tree / n_tree_subset

}





