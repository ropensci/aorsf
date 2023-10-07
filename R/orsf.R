
#' Oblique Random Survival Forest (ORSF)
#'
#' Fit an oblique random survival forest
#'
#' @srrstats {G1.4} *documented with Roxygen*
#' @srrstats {G1.1} *aorsf is an improvement of the ORSF algorithm implemented in obliqueRSF, which was an extension of Hemant Ishwaran's random survival forest.*
#' @srrstats {G1.3} *linear combinations of inputs defined.*
#' @srrstats {G1.5} *orsf() will be used in publications to benchmark performance of the aorsf package in computation speed and prediction accuracy.*
#' @srrstats {G1.6} *orsf() should be used to compare performance claims with other packages.*
#' @srrstats {G2.1} *Inputs have indication of type in parentheticals. This format is used in all exported functions.*
#' @srrstats {G5.2a} *messages produced here (e.g., with `stop()`, `warning()`, `message()`) are unique and make effort to highlight the specific data elements that cause the error*
#' @srrstats {G2.0a} *secondary documentation of arg lengths. When an input has length 1, a parenthetical gives the specific type of value it should be and uses a singular description (e.g., an integer). When inputs have length > 1, a vector description is used (e.g., integer vector)*
#' @srrstats {ML1.0} *Documentation includes a subsection that makes clear conceptual distinction between train and test data*
#' @srrstats {ML3.3} *Properties and behaviours of aorsf models are explicitly compared with objects produced by other ML software in the "Introduction to aorsf" vignette.*
#' @srrstats {ML4.0} *orsf() is a unified single-function interface to model training. orsf_train() is able to receive as input an untrained model specified by orsf() when no_fit = TRUE. Models with categorically different specifications are able to be submitted to the same model training function.*
#' @srrstats {ML5.2, ML5.2a} *The structure and functionality of trained aorsf objects is documented through vignettes. In particular, basic functionality extending from the aorsf class is explicitly described in the "Introduction to aorsf" vignette, and additional functionality is documented in the "Out-of-bag predictions and evaluation" and "Compute partial dependence with ORSF" vignettes. Each vignettes demonstrates functionality clearly with example code.*
#' @srrstats {ML5.3} *Assessment of model performance is implemented through out-of-bag error, which is finalized after a model is trained*
#' @srrstats {ML5.4} *The "Out-of-bag predictions and evaluation" vignette shows how to implement built-in or user-specified functions for this functionality.*
#' @srrstats {ML1.1} *Training data are labelled as "train".*
#' @srrstats {G2.5} *factors used as predictors can be ordered and un-ordered.*
#' @srrstats {ML4.1b} *The value of out-of-bag error can be returned for every oobag_eval_every step.*
#' @srrstats {ML4.2} *The extraction of out-of-bag error is explicitly documented with example code in the "Out-of-bag predictions and evaluation" vignette.*
#' @srrstats {ML3.5b} *Users can specify the kind of loss function to assess distance between model estimates and desired output. This is discussed in detail in the "Out-of-bag predictions and evaluation" vignette.*
#' @srrstats {ML5.4a} *Harrell's C-statistic, an internally utilized metric for model performance, is clearly and distinctly documented and cited.*
#' @srrstats {ML5.4b} *It is possible to submit custom metrics to a model assessment function, and the ability to do so is clearly documented. The "Out-of-bag predictions and evaluation" vignette provides example code.*
#' @srrstats {ML2.0, ML2.0b} *orsf() enables pre-processing steps to be defined and parametrized without fitting a model when no_fit is TRUE, returning an object with a defined class minimally intended to implement a default `print` method which summarizes the model specifications.*
#' @srrstats {ML3.0} *Model specification can be implemented prior to actual model fitting or training*
#' @srrstats {ML3.0a} *As pre-processing, model specification, and training are controlled by the orsf() function, an input parameter (no_fit) enables models to be specified yet not fitted.*
#' @srrstats {ML3.0c} *when no_fit=TRUE, orsf() will return an object that can be directly trained using orsf_train().*
#' @srrstats {ML1.6a} *Explain why missing values are not admitted.*
#' @srrstats {G1.0} *Jaeger et al describes the ORSF algorithm that aorsf is based on. Note: aorsf uses a different approach to create linear combinations of inputs for speed reasons, but orsf_control_net() allows users to make ensembles that are very similar to obliqueRSF::ORSF().*
#' @srrstats {ML1.6b} *Explicit example showing how missing values may be imputed rather than discarded.*
#' @srrstats {ML6.0} *Reference section explicitly links to aorsf-bench, which includes training and testing stages, and which clearly indicates a need for distinct training and test data sets.*
#' @srrstats {ML6.1} *clearly document how aorsf can be embedded within a typical full ML workflow.*
#' @srrstats {ML6.1a} *Embed aorsf within a full workflow using tidymodels and tidyverse*
#' @srrstats {ML5.2b} *Documentation includes examples of how to save and re-load trained model objects for their re-use.*
#' @srrstats {ML2.3} *Values associated with transformations are recorded in the object returned by orsf()*
#' @srrstats {ML1.3} *Input data are partitioned as training (in-bag) and test (out-of-bag) data within orsf_fit().*
#' @srrstats {ML4.1} *orsf_fit() retains information on model-internal parameters.*
#' @srrstats {ML4.1a} *orsf_fit() output includes all model-internal parameters, specifically the linear combination coefficients.*
#'
#' @param data a `r roxy_data_allowed()` that contains the
#'  relevant variables.
#'
#' @param formula (_formula_) The response on the left hand side should
#'   include a time variable, followed by a status variable, and may be
#'   written inside a call to [Surv][survival::Surv] (see examples).
#'   The terms on the right are names of predictor variables.
#'
#' @param control (*orsf_control*) An object returned from one of the
#'  `orsf_control` functions:
#'
#'  - [orsf_control_fast] (the default) uses a single iteration of Newton
#'    Raphson scoring to identify a linear combination of predictors.
#'
#'  - [orsf_control_cph] uses Newton Raphson scoring until a convergence
#'    criteria is met.
#'
#'  - [orsf_control_net] uses `glmnet` to identify linear combinations of
#'    predictors, similar to Jaeger (2019).
#'
#'  - [orsf_control_custom] allows the user to apply their own function
#'    to create linear combinations of predictors.
#'
#' @param weights (_numeric vector_) Optional. If given, this input should
#'   have length equal to `nrow(data)`. Values in `weights` are treated like
#'   replication weights, i.e., a value of 2 is the same thing as having 2
#'   observations in `data`, each containing a copy of the corresponding
#'   person's data.
#'
#'   *Use* `weights` *cautiously*, as `orsf` will count the number of
#'   observations and events prior to growing a node for a tree, so higher
#'   values in `weights` will lead to deeper trees.
#'
#' @param n_tree (_integer_) the number of trees to grow.
#' Default is `n_tree = 500.`
#'
#' @param n_split (_integer_) the number of cut-points assessed when splitting
#'  a node in decision trees. Default is `n_split = 5`.
#'
#' @param n_retry (_integer_) when a node can be split, but the current
#'  linear combination of inputs is unable to provide a valid split, `orsf`
#'  will try again with a new linear combination based on a different set
#'  of randomly selected predictors, up to `n_retry` times. Default is
#'  `n_retry = 3`. Set `n_retry = 0` to prevent any retries.
#'
#' @param n_thread `r roxy_n_thread_header("growing trees, computing predictions, and computing importance")`
#'
#' @param mtry (_integer_) Number of predictors randomly included as candidates
#'   for splitting a node. The default is the smallest integer greater than
#'   the square root of the number of total predictors, i.e.,
#'   `mtry = ceiling(sqrt(number of predictors))`
#'
#' @param sample_with_replacement (_logical_) If `TRUE` (the default),
#'   observations are sampled with replacement when an in-bag sample
#'   is created for a decision tree. If `FALSE`, observations are
#'   sampled without replacement and each tree will have an in-bag sample
#'   containing `sample_fraction`% of the original sample.
#'
#' @param sample_fraction (_double_) the proportion of observations that
#'   each trees' in-bag sample will contain, relative to the number of
#'   rows in `data`. Only used if `sample_with_replacement` is `FALSE`.
#'   Default value is 0.632.
#'
#' @param leaf_min_events (_integer_) minimum number of events in a
#'   leaf node. Default is `leaf_min_events = 1`
#'
#' @param leaf_min_obs (_integer_) minimum number of observations in a
#'   leaf node. Default is `leaf_min_obs = 5`.
#'
#' @param split_rule (_character_) how to assess the quality of a potential
#'   splitting rule for a node. Valid options are
#'
#'   - 'logrank' : a log-rank test statistic.
#'   - 'cstat'   : Harrell's concordance statistic.
#'
#' @param split_min_events (_integer_) minimum number of events required
#'   in a node to consider splitting it. Default is `split_min_events = 5`
#'
#' @param split_min_obs (_integer_) minimum number of observations required
#'   in a node to consider splitting it. Default is `split_min_obs = 10`.
#'
#' @param split_min_stat (double) minimum test statistic required to split
#'   a node. Default is 3.841459 if `split_rule = 'logrank'` and 0.50 if
#'   `split_rule = 'cstat'`. If no splits are found with a statistic
#'   exceeding `split_min_stat`, the given node either becomes a leaf or
#'   a retry occurs (up to `n_retry` retries).
#'
#' @param oobag_pred_type (_character_) The type of out-of-bag predictions
#'   to compute while fitting the ensemble. Valid options are
#'
#'   - 'none' : don't compute out-of-bag predictions
#'   - 'risk' : probability of event occurring at or before `oobag_pred_horizon`.
#'   - 'surv' : 1 - risk.
#'   - 'chf'  : cumulative hazard function at `oobag_pred_horizon`.
#'   - 'mort' : mortality, i.e., the number of events expected if all
#'              observations in the training data were identical to a
#'              given observation.
#'
#' @param oobag_pred_horizon (_numeric_) A numeric value indicating what time
#'   should be used for out-of-bag predictions. Default is the median
#'   of the observed times, i.e., `oobag_pred_horizon = median(time)`.
#'
#' @param oobag_eval_every (_integer_) The out-of-bag performance of the
#'   ensemble will be checked every `oobag_eval_every` trees. So, if
#'   `oobag_eval_every = 10`, then out-of-bag performance is checked
#'   after growing the 10th tree, the 20th tree, and so on. Default
#'   is `oobag_eval_every = n_tree`.
#'
#' @param oobag_fun `r roxy_oobag_fun_header()` every `oobag_eval_every`
#'  trees. `r roxy_oobag_fun_default()` `r roxy_oobag_fun_user()`
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
#' @param group_factors (_logical_) Only relevant if variable importance is
#'   being estimated. `r roxy_group_factors()`
#'
#' @param tree_seeds (_integer vector_) Optional. if specified, random seeds
#'   will be set using the values in `tree_seeds[i]`  before growing tree `i`.
#'   Two forests grown with the same number of trees and the same seeds will
#'   have the exact same out-of-bag samples, making out-of-bag error
#'   estimates of the forests more comparable. If `NULL` (the default),
#'   no seeds are set during the training process.
#'
#' @param attach_data (_logical_) if `TRUE`, a copy of the training
#'   data will be attached to the output. This is helpful if you
#'   plan on using functions like [orsf_pd_oob] or [orsf_summarize_uni]
#'   to interpret the forest using its training data. Default is `TRUE`.
#'
#' @param no_fit (_logical_) if `TRUE`, model fitting steps are defined and
#'   saved, but training is not initiated. The object returned can be
#'   directly submitted to `orsf_train()` so long as `attach_data` is `TRUE`.
#'
#' @param na_action `r roxy_na_action_header("data")`
#'
#'   - `r roxy_na_action_fail("data")`
#'   - `r roxy_na_action_omit("data")`
#'   - `r roxy_na_action_impute_meanmode("data")`. Note that is this
#'     option is selected and `attach_data` is `TRUE`, the data attached
#'     to the output will be the imputed version of `data`.
#'
#' @param verbose_progress (_logical_) if `TRUE`, progress messages are
#'   printed in the console. If `FALSE` (the default), nothing is printed.
#'
#' @param ... `r roxy_dots()`
#'
#' @param object an untrained 'aorsf' object, created by setting
#'   `no_fit = TRUE` in `orsf()`.
#'
#' @srrstats {ML5.0} *The result of applying orsf training processes results in a single model object.*
#'
#' @return an accelerated oblique RSF object (`aorsf`)
#'
#' @details
#'
#' This function is based on and similar to the `ORSF` function
#'   in the `obliqueRSF` R package. The primary difference is that this
#'   function runs much faster. The speed increase is attributable to better
#'   management of memory (i.e., no unnecessary copies of inputs) and using
#'   a Newton Raphson scoring algorithm to identify linear combinations of
#'   inputs rather than performing penalized regression using routines in
#'   `glmnet`.The modified Newton Raphson scoring algorithm that this
#'   function applies is an adaptation of the C++ routine developed by
#'   Terry M. Therneau that fits Cox proportional hazards models
#'   (see [survival::coxph()] and more specifically [survival::coxph.fit()]).
#'
#'
#' @section Details on inputs:
#'
#' **formula**:
#'
#' - The response in `formula` can be a survival
#'   object as returned by the [Surv][survival::Surv] function,
#'   but can also just be the time and status variables. I.e.,
#'   `Surv(time, status) ~ .` works just like `time + status ~ .`
#'
#' - A `.` symbol on the right hand side is short-hand for using all
#'   variables in `data` (omitting those on the left hand side of
#'   `formula`) as predictors.
#'
#' - The order of variables in the left hand side matters. i.e.,
#'   writing `status + time ~ .` will make `orsf` assume your
#'   `status` variable is actually the `time` variable.
#'
#' - The response variable can be a survival object stored in `data`.
#'   For example, y ~ . is a valid formula if `data$y` inherits
#'   from the `Surv` class.
#'
#' - Although you can fit an oblique random survival forest with 1 predictor
#'   variable, your formula should have at least 2 predictors. The reason for
#'   this recommendation is that a linear combination of predictors is trivial
#'   if there is only one predictor.
#'
#'
#' **mtry**:
#'
#' The `mtry` parameter may be temporarily reduced to ensure there
#'   are at least 2 events per predictor variable. This occurs when using
#'   [orsf_control_cph] because coefficients in the Newton Raphson scoring
#'   algorithm may become unstable when the number of covariates is
#'   greater than or equal to the number of events. This reduction does not
#'   occur when using [orsf_control_net].
#'
#' **oobag_fun**:
#'
#' If `oobag_fun` is specified, it will be used in to compute negation
#'  importance or permutation importance, but it will not have any role
#'  for ANOVA importance.
#'
#' **n_thread**:
#'
#' If an R function must be called from C++ (i.e., user-supplied function to
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
#' @section Missing data:
#'
#' Data passed to aorsf functions are not allowed to have missing values.
#'   A user should impute missing values using an R package with that purpose,
#'   such as `recipes` or `mlr3pipelines`.
#'
#' @references
#'
#' `r roxy_cite_harrell_1982()`
#'
#' `r roxy_cite_breiman_2001()`
#'
#' `r roxy_cite_ishwaran_2008()`
#'
#' `r roxy_cite_jaeger_2019()`
#'
#' `r roxy_cite_jaeger_2022()`
#'
#' @export
#'
#' @includeRmd Rmd/orsf_examples.Rmd
#'
#'

orsf <- function(data,
                 formula,
                 control = orsf_control_fast(),
                 weights = NULL,
                 n_tree = 500,
                 n_split = 5,
                 n_retry = 3,
                 n_thread = 1,
                 mtry = NULL,
                 sample_with_replacement = TRUE,
                 sample_fraction = 0.632,
                 leaf_min_events = 1,
                 leaf_min_obs = 5,
                 split_rule = 'logrank',
                 split_min_events = 5,
                 split_min_obs = 10,
                 split_min_stat = switch(split_rule,
                                         "logrank" = 3.841459,
                                         "cstat" = 0.50),
                 oobag_pred_type = 'surv',
                 oobag_pred_horizon = NULL,
                 oobag_eval_every = n_tree,
                 oobag_fun = NULL,
                 importance = 'anova',
                 group_factors = TRUE,
                 tree_seeds = NULL,
                 attach_data = TRUE,
                 no_fit = FALSE,
                 na_action = 'fail',
                 verbose_progress = FALSE,
                 ...){

 #' @srrstats {G2.8} *As part of initial pre-processing, run checks on inputs to ensure that all other sub-functions receive inputs of a single defined class or type.*

 check_dots(list(...), .f = orsf)

 if(!is.data.frame(data))
  data <- orsf_data_prep(data)

 check_orsf_inputs(
  data = data,
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
  oobag_pred_type = oobag_pred_type,
  oobag_pred_horizon = oobag_pred_horizon,
  oobag_eval_every = oobag_eval_every,
  importance = importance,
  tree_seeds = tree_seeds,
  attach_data = attach_data,
  verbose_progress = verbose_progress
 )

 #TODO: more polish
 if(split_rule == "cstat" && split_min_stat >= 1){
  stop("If split_rule is 'cstat', split_min_stat must be < 1",
       call. = FALSE)
 }

 oobag_pred <- oobag_pred_type != 'none'

 if(sample_fraction == 1 && oobag_pred){
  stop(
   "cannot compute out-of-bag predictions if no samples are out-of-bag.",
   "To resolve this, set sample_fraction < 1 or oobag_pred_type = 'none'.",
   call. = FALSE
  )
 }

 orsf_type <- attr(control, 'type')

 switch(
  orsf_type,

  'fast' = {

   control_net <- orsf_control_net()
   control_cph <- control
   f_beta      <- function(x) x

  },

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
   control_cph <- orsf_control_fast(do_scale = FALSE)
   f_beta      <- penalized_cph
  },

  "custom" = {

   control_net <- orsf_control_net()
   control_cph <- orsf_control_fast(do_scale = FALSE)
   f_beta      <- control$beta_fun

  }

 )

 # TODO: drop this; importance is computed with mort
 if(importance %in% c("permute", "negate") && !oobag_pred){
  # oobag_pred <- TRUE # Should I add a warning?
  oobag_pred_type <- 'surv'
 }

 if(is.null(oobag_fun)){

  f_oobag_eval <- function(x) x
  type_oobag_eval <- if(oobag_pred) 'cstat' else 'none'

 } else {

  check_oobag_fun(oobag_fun)
  f_oobag_eval <- oobag_fun
  type_oobag_eval <- 'user'

  if(oobag_pred_type == 'leaf'){
   stop("a user-supplied oobag function cannot be",
        "applied when oobag_pred_type = 'leaf'",
        call. = FALSE)
  }

 }

 # can't evaluate the oobag predictions if they aren't aggregated
 if(oobag_pred_type == 'leaf') type_oobag_eval <- 'none'


 cph_method <- control_cph$cph_method
 cph_eps <- control_cph$cph_eps
 cph_iter_max <- control_cph$cph_iter_max
 cph_do_scale <- control_cph$cph_do_scale
 net_alpha <- control_net$net_alpha
 net_df_target <- control_net$net_df_target


 formula_terms <- suppressWarnings(stats::terms(formula, data=data))

 if(attr(formula_terms, 'response') == 0)
  stop("formula must have a response", call. = FALSE)

 names_y_data <- all.vars(formula[[2]])

 if(length(names_y_data) == 1){
  # this is fine if the response is a Surv object,
  if(!inherits(data[[names_y_data]], 'Surv')){
   # otherwise it will be a problem
   stop("formula must have two variables (time & status) as the response",
        call. = FALSE)
  }

 }

 if(length(names_y_data) > 2){
  stop("formula must have two variables (time & status) as the response",
       call. = FALSE)
 }

 types_y_data <- vector(mode = 'character',
                        length = length(names_y_data))

 for(i in seq_along(types_y_data)){
  types_y_data[i] <- class(data[[ names_y_data[i] ]])[1]
 }

 unit_y_names <- names_y_data[types_y_data == 'units']

 ui_y <- unit_info(data = data, .names = unit_y_names)

 names_x_data <- attr(formula_terms, 'term.labels')

 names_not_found <- setdiff(c(names_y_data, names_x_data), names(data))

 if(!is_empty(names_not_found)){
  msg <- paste0(
   "variables in formula were not found in data: ",
   paste_collapse(names_not_found, last = ' and ')
  )
  stop(msg, call. = FALSE)
 }

 #' @srrstats {G2.7} *aorsf accepts as input numeric and categorical predictor variables, including those with unit class. I do not think it is necessary to incorporate any other type of input, since it is relatively straightforward to convert data into a numeric or categorical format.*

 #' @srrstats {G2.11} *I cannot write code for every possible vector class to ensure that the vector data will be safely coerced into a a valid class and the attributes will be stored in the orsf_out object. It is much easier and safer for the user to convert a few columns to numeric or factor than it is for me to attempt writing code that will safely coerce every type of vector. That being said, I do find units columns to be helpful and I've written some code to make them an allowable class in input data.*

 types_x_data <- check_var_types(data,
                                 names_x_data,
                                 valid_types = c('numeric',
                                                 'integer',
                                                 'units',
                                                 'factor',
                                                 'ordered'))

 #' @srrstats {G2.6} *ensure that one-dimensional inputs are appropriately pre-processed. aorsf does not deal with missing data as many other R packages are very good at dealing with it.*

 #' @srrstats {G2.13} *check for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*

 #' @srrstats {G2.15} *Never pass data with potential missing values to any base routines.*

 #' @srrstats {G2.16} *Throw hard errors if undefined values are detected.*

 for(i in c(names_y_data, names_x_data)){

  if(any(is.infinite(data[[i]]))){
   stop("Please remove infinite values from ", i, ".",
        call. = FALSE)
  }

  if(any(collapse::allNA(data[[i]]))){
   stop("column ", i, " has no observed values",
        call. = FALSE)
  }

  # nan values trigger is.na(), so this probably isnt needed.
  # if(any(is.nan(data[[i]]))){
  #  stop("Please remove NaN values from ", i, ".",
  #       call. = FALSE)
  # }

 }

 fctr_check(data, names_x_data)

 fctr_id_check(data, names_x_data)

 fi <- fctr_info(data, names_x_data)

 unit_x_names <- names_x_data[types_x_data == 'units']

 ui_x <- unit_info(data = data, .names = unit_x_names)

 names_x_numeric <- grep(pattern = "^integer$|^numeric$|^units$",
                         x = types_x_data)

 means <- standard_deviations<- modes <- numeric_bounds <- NULL

 numeric_cols <- names_x_data[names_x_numeric]
 nominal_cols <- fi$cols

 if(!is_empty(nominal_cols)){

  modes <- vapply(
   select_cols(data, nominal_cols),
   collapse::fmode,
   FUN.VALUE = integer(1),
   w = weights
  )

 }

 if(!is_empty(numeric_cols)){

  numeric_data <- select_cols(data, numeric_cols)

  numeric_bounds <- matrix(
   data = c(
    collapse::fnth(numeric_data, 0.1),
    collapse::fnth(numeric_data, 0.25),
    collapse::fnth(numeric_data, 0.5),
    collapse::fnth(numeric_data, 0.75),
    collapse::fnth(numeric_data, 0.9)
   ),
   nrow =5,
   byrow = TRUE,
   dimnames = list(c('10%', '25%', '50%', '75%', '90%'),
                   names(numeric_data))
  )

  means <- collapse::fmean(numeric_data, w = weights)

  standard_deviations <- collapse::fsd(numeric_data, w = weights)

 }

 if(any(is.na(select_cols(data, names_y_data))))
  stop("Please remove missing values from the outcome variable(s)",
       call. = FALSE)

 if(any(is.na(select_cols(data, names_x_data)))){

  switch(
   na_action,

   'fail' = {
    stop("Please remove missing values from data, or impute them.",
         call. = FALSE)
   },

   'omit' = {
    # data <- collapse::na_omit(data, cols = names_x_data)
    cc <- stats::complete.cases(data[, names_x_data])
    data <- data[cc, ]
   },

   'impute_meanmode' = {

    data <- data_impute(data,
                        cols = names_x_data,
                        values = c(as.list(means),
                                   as.list(modes)))
   }

  )

 }

 y <- prep_y(data, names_y_data)
 x <- prep_x(data, fi, names_x_data, means, standard_deviations)

 if(is.null(mtry)) mtry <- ceiling(sqrt(ncol(x)))

 if(is.null(net_df_target)) net_df_target <- mtry

 # warn instead?
 if(net_df_target > mtry)
  stop("net_df_target = ", net_df_target,
       " must be <= mtry, which is ", mtry,
       call. = FALSE)

 n_events <- collapse::fsum(y[, 2])

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

 if(!is.null(oobag_pred_horizon)){

  if(any(oobag_pred_horizon <= 0))

   stop("Out of bag prediction horizon (oobag_pred_horizon) must be > 0",
        call. = FALSE)

 } else {

  # use training data to provide sensible default
  oobag_pred_horizon <- stats::median(y[, 1])

 }

 sorted <-
  collapse::radixorder(y[, 1],  # order this way for risk sets
                       -y[, 2]) # order this way for oob C-statistic.

 if(is.null(weights)) weights <- rep(1, nrow(x))

 x_sort <- x[sorted, , drop = FALSE]
 y_sort <- y[sorted, , drop = FALSE]
 w_sort <- weights[sorted]

 if(length(tree_seeds) == 1 && n_tree > 1){
  set.seed(tree_seeds)
  tree_seeds <- sample(x = n_tree*10, size = n_tree, replace = FALSE)
 } else if(is.null(tree_seeds)){
  tree_seeds <- sample(x = n_tree*10, size = n_tree, replace = FALSE)
 }


 vi_max_pvalue = 0.01
 tree_type_R = 3

 orsf_out <- orsf_cpp(x = x_sort,
                      y = y_sort,
                      w = w_sort,
                      tree_type_R = tree_type_R,
                      tree_seeds = as.integer(tree_seeds),
                      loaded_forest = list(),
                      n_tree = n_tree,
                      mtry = mtry,
                      sample_with_replacement = sample_with_replacement,
                      sample_fraction = sample_fraction,
                      vi_type_R = switch(importance,
                                         "none" = 0,
                                         "negate" = 1,
                                         "permute" = 2,
                                         "anova" = 3),
                      vi_max_pvalue = vi_max_pvalue,
                      lincomb_R_function = f_beta,
                      oobag_R_function = f_oobag_eval,
                      leaf_min_events = leaf_min_events,
                      leaf_min_obs = leaf_min_obs,
                      split_rule_R = switch(split_rule,
                                            "logrank" = 1,
                                            "cstat" = 2),
                      split_min_events = split_min_events,
                      split_min_obs = split_min_obs,
                      split_min_stat = split_min_stat,
                      split_max_cuts = n_split,
                      split_max_retry = n_retry,
                      lincomb_type_R = switch(orsf_type,
                                              'fast' = 1,
                                              'cph' = 1,
                                              'random' = 2,
                                              'net' = 3,
                                              'custom' = 4),
                      lincomb_eps = cph_eps,
                      lincomb_iter_max = cph_iter_max,
                      lincomb_scale = cph_do_scale,
                      lincomb_alpha = net_alpha,
                      lincomb_df_target = net_df_target,
                      lincomb_ties_method = switch(tolower(cph_method),
                                                   'breslow' = 0,
                                                   'efron'   = 1),
                      pred_type_R = switch(oobag_pred_type,
                                           "none" = 0,
                                           "risk" = 1,
                                           "surv" = 2,
                                           "chf"  = 3,
                                           "mort" = 4,
                                           "leaf" = 8),
                      pred_mode = FALSE,
                      pred_aggregate = oobag_pred_type != 'leaf',
                      pred_horizon = oobag_pred_horizon,
                      oobag = oobag_pred,
                      oobag_eval_type_R = switch(type_oobag_eval,
                                                 'none' = 0,
                                                 'cstat' = 1,
                                                 'user' = 2),
                      oobag_eval_every = oobag_eval_every,
                      pd_type_R = 0,
                      pd_x_vals = list(matrix(0, ncol=1, nrow=1)),
                      pd_x_cols = list(matrix(1L, ncol=1, nrow=1)),
                      pd_probs = c(0),
                      n_thread = n_thread,
                      write_forest = TRUE,
                      run_forest = !no_fit,
                      verbosity = as.integer(verbose_progress))

 # if someone says no_fit and also says don't attach the data,
 # give them a warning but also do the right thing for them.
 orsf_out$data <- if(attach_data) data else NULL

 if(importance != 'none' && !no_fit){
  rownames(orsf_out$importance) <- colnames(x)
  orsf_out$importance <-
   rev(orsf_out$importance[order(orsf_out$importance), , drop=TRUE])
 }

 if(oobag_pred){

  # makes labels for oobag evaluation type
  orsf_out$eval_oobag$stat_type <-
   switch(EXPR = as.character(orsf_out$eval_oobag$stat_type),
          "0" = "None",
          "1" = "Harrell's C-statistic",
          "2" = "User-specified function")

  if(!no_fit){

  # put the oob predictions into the same order as the training data.
  unsorted <- collapse::radixorder(sorted)

  if(oobag_pred_type == 'leaf'){
   all_rows <- seq(nrow(data))
   for(i in seq(n_tree)){
    rows_inbag <- setdiff(all_rows, orsf_out$forest$rows_oobag[[i]]+1)
    orsf_out$pred_oobag[rows_inbag, i] <- NA
   }
  }

  orsf_out$pred_oobag <- orsf_out$pred_oobag[unsorted, , drop = FALSE]
  orsf_out$pred_oobag[is.nan(orsf_out$pred_oobag)] <- NA_real_

  # mortality predictions should always be 1 column
  # b/c they do not depend on the prediction horizon
  if(oobag_pred_type == 'mort'){
   orsf_out$pred_oobag <-
    orsf_out$pred_oobag[, 1L, drop = FALSE]

   orsf_out$eval_oobag$stat_values <-
    orsf_out$eval_oobag$stat_values[, 1, drop = FALSE]
  }

  }

 }

 orsf_out$pred_horizon <- oobag_pred_horizon

 n_leaves_mean <- compute_mean_leaves(orsf_out$forest)

 attr(orsf_out, 'control')             <- control
 attr(orsf_out, 'mtry')                <- mtry
 attr(orsf_out, 'n_obs')               <- nrow(y_sort)
 attr(orsf_out, 'n_tree')              <- n_tree
 attr(orsf_out, 'names_y')             <- names_y_data
 attr(orsf_out, "names_x")             <- names_x_data
 attr(orsf_out, "names_x_ref")         <- colnames(x)
 attr(orsf_out, "types_x")             <- types_x_data
 attr(orsf_out, 'n_events')            <- n_events
 attr(orsf_out, 'max_time')            <- y_sort[nrow(y_sort), 1]
 attr(orsf_out, 'event_times')         <- unique(y_sort[ y_sort[,2]==1, 1])
 attr(orsf_out, "unit_info")           <- c(ui_y, ui_x)
 attr(orsf_out, "fctr_info")           <- fi
 attr(orsf_out, 'n_leaves_mean')       <- n_leaves_mean
 attr(orsf_out, 'n_split')             <- n_split
 attr(orsf_out, 'leaf_min_events')     <- leaf_min_events
 attr(orsf_out, 'leaf_min_obs')        <- leaf_min_obs
 attr(orsf_out, 'split_min_events')    <- split_min_events
 attr(orsf_out, 'split_min_obs')       <- split_min_obs
 attr(orsf_out, 'split_min_stat')      <- split_min_stat
 attr(orsf_out, 'na_action')           <- na_action
 attr(orsf_out, 'cph_method')          <- cph_method
 attr(orsf_out, 'cph_eps')             <- cph_eps
 attr(orsf_out, 'cph_iter_max')        <- cph_iter_max
 attr(orsf_out, 'cph_do_scale')        <- cph_do_scale
 attr(orsf_out, 'net_alpha')           <- net_alpha
 attr(orsf_out, 'net_df_target')       <- net_df_target
 attr(orsf_out, 'numeric_bounds')      <- numeric_bounds
 attr(orsf_out, 'means')               <- means
 attr(orsf_out, 'modes')               <- modes
 attr(orsf_out, 'standard_deviations') <- standard_deviations
 attr(orsf_out, 'trained')             <- !no_fit
 attr(orsf_out, 'n_retry')             <- n_retry
 attr(orsf_out, 'orsf_type')           <- orsf_type
 attr(orsf_out, 'f_beta')              <- f_beta
 attr(orsf_out, 'f_oobag_eval')        <- f_oobag_eval
 attr(orsf_out, 'type_oobag_eval')     <- type_oobag_eval
 attr(orsf_out, 'oobag_pred')          <- oobag_pred
 attr(orsf_out, 'oobag_fun')           <- oobag_fun
 attr(orsf_out, 'oobag_pred_type')     <- oobag_pred_type
 attr(orsf_out, 'oobag_eval_every')    <- oobag_eval_every
 attr(orsf_out, 'oobag_pred_horizon')  <- oobag_pred_horizon
 attr(orsf_out, 'importance')          <- importance
 attr(orsf_out, 'importance_values')   <- orsf_out$importance
 attr(orsf_out, 'group_factors')       <- group_factors
 attr(orsf_out, 'weights_user')        <- weights
 attr(orsf_out, 'verbose_progress')    <- verbose_progress
 attr(orsf_out, 'vi_max_pvalue')       <- vi_max_pvalue
 attr(orsf_out, 'split_rule')          <- split_rule
 attr(orsf_out, 'n_thread')            <- n_thread
 attr(orsf_out, 'tree_type')           <- tree_type_R
 attr(orsf_out, 'tree_seeds')          <- tree_seeds
 attr(orsf_out, 'sample_with_replacement') <- sample_with_replacement
 attr(orsf_out, 'sample_fraction')         <- sample_fraction

 #' @srrstats {ML5.0a} *orsf output has its own class*
 class(orsf_out) <- "orsf_fit"

 if(importance != 'none' && !no_fit){

  # temporarily attach data
  if(!attach_data) orsf_out$data <- data

  orsf_out$importance <- orsf_vi(orsf_out,
                                 importance = importance,
                                 group_factors = group_factors)

  # Can drop it now
  if(!attach_data) orsf_out$data <- NULL

 }

 orsf_out

}


orsf_data_prep <- function(data, ...){
 UseMethod('orsf_data_prep')
}

orsf_data_prep.list <- function(data, ...){

 lengths <- vapply(data, length, integer(1))

 if(! all(lengths == lengths[1])){

  length_tbl <- table(lengths)
  length_mode <- as.numeric(names(length_tbl)[which.max(length_tbl)])

  mismatch <- lengths[names(which(lengths != length_mode))]

  mismatch <-
   paste(" -", names(mismatch),
         'has length', mismatch,
         collapse = '\n')

  mismatch <-
   paste(mismatch,
         '\n - all other variables have length ', length_mode,
         sep = '')

  stop("unable to cast data (a list) into a data.frame.\n",
       mismatch, call. = FALSE)

 }

 data <-
  tryCatch(as.data.frame(data), error = function(e) e$message)

 if(!is.data.frame(data)){
  stop("Could not coerce data (a list) into a data.frame object.\n",
       "Running as.data.frame(data) ",
       "produced this error message:\n\"", data, "\"",
       call. = FALSE)
 }

 data

}

orsf_data_prep.recipe <- function(data, ...){

 getElement(data, 'template')

}

#' @srrstats {ML2.0a} *objects returned from orsf() with no_fit = TRUE can be directly submitted to orsf_train to train the model specification.*
#'
#' @rdname orsf
#' @export
orsf_train <- function(object){
 orsf_train_(object)
}



#' Estimate training time
#'
#' @srrstats {ML4.5} *include a function to estimate likely time to train a specified model*
#'
#' @param object an untrained `aorsf` object
#'
#' @param n_tree_subset (_integer_)  how many trees should be fit in order
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
#'                n_tree = 500, no_fit = TRUE)
#'
#' # grow 50 trees to approximate the time it will take to grow 500 trees
#' time_estimated <- orsf_time_to_train(object, n_tree_subset = 50)
#'
#' print(time_estimated)
#'
#' # let's see how close the approximation was
#' time_true_start <- Sys.time()
#' fit <- orsf_train(object)
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

 time_preproc_start <- Sys.time()

 y  <- prep_y_from_orsf(object)
 x  <- prep_x_from_orsf(object)

 sorted <- order(y[, 1],  # order this way for risk sets
                 -y[, 2]) # order this way for oob C-statistic.

 time_preproc_stop <- Sys.time()

 output <- orsf_train_(object,
                       n_tree = n_tree_subset,
                       x = x,
                       y = y,
                       sorted = sorted)

 time_train_stop <- Sys.time()

 mult_by <- get_n_tree(object) / n_tree_subset

 difftime(time_preproc_stop, time_preproc_start) +
  difftime(time_train_stop, time_preproc_stop) * mult_by

}


#' internal training function
#'
#' the purpose of this function is to have more control over fit parameters
#'   so that I can do tasks other than train the specific aorsf specification,
#'   e.g., estimate how much time it will take to train it.
#'
#' @param object see orsf()
#' @param n_tree default is NULL. If this input is specified, then it will
#'   be used instead of the n_tree attribute in object.
#' @param x default is NULL. If this input is specified, then it will
#'   be used instead of creating an x matrix by normal means.
#' @param y default is NULL. If this input is specified, then it will
#'   be used instead of creating a y matrix by normal means.
#' @param sorted default is NULL. If this input is specified, then it will
#'   be used instead of sorting the y matrix.
#'
#' @return a trained aorsf model
#'
#' @noRd
#'
orsf_train_ <- function(object,
                        n_tree = NULL,
                        x = NULL,
                        y = NULL,
                        sorted  = NULL){

 if(is_trained(object)){
  stop("object has already been trained", call. = FALSE)
 }

 if(is.null(object$data)){
  stop("object must have training data attached.",
       " Set attach_data = TRUE in orsf()",
       call. = FALSE)
 }

 if(is.null(y)){
  y <- prep_y_from_orsf(object)
 }

 if(is.null(x)){
  x <- prep_x_from_orsf(object)
 }

 if(is.null(n_tree)){
  n_tree <- get_n_tree(object)
 }

 if(is.null(sorted)){
  sorted <- collapse::radixorder(y[, 1], -y[, 2])
 }

 weights <- get_weights_user(object)

 x_sort <- x[sorted, ]
 y_sort <- y[sorted, ]
 w_sort <- weights[sorted]

 oobag_eval_every <- min(n_tree, get_oobag_eval_every(object))

 orsf_out <- orsf_cpp(x = x_sort,
                      y = y_sort,
                      w = w_sort,
                      tree_type_R = 3,
                      tree_seeds = get_tree_seeds(object),
                      loaded_forest = list(),
                      n_tree = n_tree,
                      mtry = get_mtry(object),
                      sample_with_replacement = get_sample_with_replacement(object),
                      sample_fraction = get_sample_fraction(object),
                      vi_type_R = switch(get_importance(object),
                                         "none" = 0,
                                         "negate" = 1,
                                         "permute" = 2,
                                         "anova" = 3),
                      vi_max_pvalue = get_vi_max_pvalue(object),
                      lincomb_R_function = get_f_beta(object),
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
                      lincomb_type_R = switch(get_orsf_type(object),
                                              'fast' = 1,
                                              'cph' = 1,
                                              'random' = 2,
                                              'net' = 3,
                                              'custom' = 4),
                      lincomb_eps = get_cph_eps(object),
                      lincomb_iter_max = get_cph_iter_max(object),
                      lincomb_scale = get_cph_do_scale(object),
                      lincomb_alpha = get_net_alpha(object),
                      lincomb_df_target = get_net_df_target(object),
                      lincomb_ties_method = switch(
                       tolower(get_cph_method(object)),
                       'breslow' = 0,
                       'efron'   = 1
                      ),
                      pred_type_R = switch(get_oobag_pred_type(object),
                                           "none" = 0,
                                           "risk" = 1,
                                           "surv" = 2,
                                           "chf"  = 3,
                                           "mort" = 4),
                      pred_mode = FALSE,
                      pred_aggregate = get_oobag_pred_type(object) != 'leaf',
                      pred_horizon = get_oobag_pred_horizon(object),
                      oobag = get_oobag_pred(object),
                      oobag_eval_type_R = switch(get_type_oobag_eval(object),
                                                 'none' = 0,
                                                 'cstat' = 1,
                                                 'user' = 2),
                      oobag_eval_every = oobag_eval_every,
                      pd_type_R = 0,
                      pd_x_vals = list(matrix(0, ncol=1, nrow=1)),
                      pd_x_cols = list(matrix(1L, ncol=1, nrow=1)),
                      pd_probs = c(0),
                      n_thread = get_n_thread(object),
                      write_forest = TRUE,
                      run_forest = TRUE,
                      verbosity = get_verbose_progress(object))


 object$pred_oobag   <- orsf_out$pred_oobag
 object$eval_oobag   <- orsf_out$eval_oobag
 object$forest       <- orsf_out$forest
 object$importance   <- orsf_out$importance
 object$pred_horizon <- get_oobag_pred_horizon(object)

 if(get_importance(object) != 'none'){

  rownames(object$importance) <- colnames(x)

  object$importance <-
   rev(object$importance[order(object$importance), , drop=TRUE])

  attr(object, 'importance_values') <- object$importance

  object$importance <- orsf_vi(object,
                               importance = get_importance(object),
                               group_factors = get_group_factors(object))


 }

 if(get_oobag_pred(object)){

  # put the oob predictions into the same order as the training data.
  # TODO: this can be faster; see predict unsorting
  unsorted <- vector(mode = 'integer', length = length(sorted))
  for(i in seq_along(unsorted)) unsorted[ sorted[i] ] <- i

  # clear labels for oobag evaluation type

  object$eval_oobag$stat_type <-
   switch(EXPR = as.character(object$eval_oobag$stat_type),
          "0" = "None",
          "1" = "Harrell's C-statistic",
          "2" = "User-specified function")

  object$pred_oobag <- object$pred_oobag[unsorted, , drop = FALSE]

  # mortality predictions should always be 1 column
  # b/c they do not depend on the prediction horizon
  if(get_oobag_pred_type(object) == 'mort'){
   object$pred_oobag <-
    object$pred_oobag[, 1L, drop = FALSE]

   object$eval_oobag$stat_values <-
    object$eval_oobag$stat_values[, 1, drop = FALSE]
  }

 }

 attr(object, "n_leaves_mean") <- compute_mean_leaves(orsf_out$forest)

 attr(object, 'trained') <- TRUE

 object

}


