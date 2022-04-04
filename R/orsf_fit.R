

#' Oblique Random Survival Forest (ORSF)
#'
#' @srrstats {G1.4} *documented with Roxygen*
#'
#' @srrstats {G1.1} *aorsf is an improvement of the ORSF algorithm implemented in obliqueRSF, which was an extension of Hemant Ishwaran's random survival forest.*
#'
#' @srrstats {G1.3} *linear combinations of inputs defined.*
#'
#' @srrstats {G1.5} *orsf() will be used in publications to benchmark performance of the aorsf package in computation speed and prediction accuracy.*
#'
#' @srrstats {G1.6} *orsf() should be used to compare performance claims with other packages.*
#'
#'
#' @srrstats {G2.1} *Inputs have indicatigiton of type in parentheticals. This format is used in all exported functions.*
#'
#' @srrstats {G5.2a} *messages produced here (e.g., with `stop()`, `warning()`, `message()`) are unique and make effort to highlight the specific data elements that cause the error*
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
#' @srrstats {ML1.1} *Training data are labelled as "train".*
#'
#' @param data_train (_data.frame_) that will be used to grow the forest.
#'
#' @srrstats {G2.5} factors used as predictors can be ordered and un-ordered.
#'
#' @param formula (_formula_) a formula object, with the response on the left
#'   of a `~` operator, and the terms on the right (see details). Variables
#'   on the right hand size of the `~` can be numeric, integer, or factor
#'   variables. Factors may be ordered or unordered.
#'
#' @param control An `aorsf_control` object, created with [orsf_control_net]
#'  or [orsf_control_cph]. Default is `control = orsf_control_cph()`.
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
#'  of randomly selected predictors, up to `n_retry` pred_horizon. When
#'  `n_retry = 0` the retry mechanic is not applied.
#'  Default is `n_retry = 0`.
#'
#' @param mtry (_integer_) Number of variables randomly selected as candidates
#'   for splitting a node. The default is the smallest integer greater than
#'   the square root of the number of features, i.e.,
#'   `mtry = ceiling(sqrt(number of predictors))`
#'
#' @param leaf_min_events (_integer_) minimum number of events in a
#'   leaf node. Default is `leaf_min_events = 1`
#'
#' @param leaf_min_obs (_integer_) minimum number of observations in a
#'   leaf node. Default is `leaf_min_obs = 5`
#'
#' @param split_min_events (_integer_) minimum number of events required
#'   to split a node. Default is `split_min_events = 5`
#'
#' @param split_min_obs (_integer_) minimum number of observations required
#'   to split a node. Default is `split_min_obs = 10`.
#'
#' @param oobag_pred (_logical_) if `TRUE` out-of-bag predictions are returned
#'   in the `aorsf` object. Default is `TRUE`.
#'
#' @param oobag_time (_numeric_) A numeric value indicating what time
#'   should be used for out-of-bag predictions. Default is the median
#'   of the observed pred_horizon, i.e., `oobag_time = median(time)`.
#'
#' @param oobag_eval_every (_integer_) The out-of-bag performance of the
#'   ensemble will be checked every `oobag_eval_every` trees. So, if
#'   `oobag_eval_every = 10`, then out-of-bag performance is checked
#'   after growing the 10th tree, the 20th tree, and so on. Default
#'   is `oobag_eval_every = n_tree`, so that out-of-bag performance is
#'   assessed once after growing all the trees.
#'
#' @param oobag_fun (_function_) When `oobag_fun` = `NULL` (the default),
#'   out-of-bag predictions are evaluated using Harrell's C-statistic.
#'   If a value for `oobag_fun` is provided, it will be used in place of
#'   Harrell's C-statistic to evaluate out-of-bag predictions. The function
#'   must have two inputs: `y_mat` and `s_vec`. The input `y_mat` is
#'   presumed to be a matrix with two columns named `time` (first column)
#'   and `status` (second column). The input `s_vec` is presumed to be a
#'   numeric vector containing predicted survival probabilities for `y_mat`.
#'
#' @param importance (_logical_) if `TRUE`, variable importance will be
#'   computed using _negation_ importance. With negation importance,
#'   all coefficients for a given variable are multiplied by -1 and
#'   then the out-of-bag error for the forest is re-computed. The greater
#'   the degradation of the forest's error, the more important the variable.
#'   Default is `FALSE`. Note that if `oobag_fun` is specified above, it
#'   will be used in the computation of negation importance.
#'
#' @param attach_data (_logical_) if `TRUE`, a copy of the training
#'   data will be attached to the output. This is helpful if you
#'   plan on using functions like [orsf_pd_summary] to interpret the fitted
#'   forest using its training data. Default is `TRUE`.
#'
#' @srrstats {ML2.0, ML2.0b} *orsf() enables pre-processing steps to be defined and parametrized without fitting a model when no_fit is TRUE, returning an object with a defined class minimally intended to implement a default `print` method which summarizes the model specifications.*
#'
#' @srrstats {ML3.0} *Model specification can be implemented prior to actual model fitting or training*
#' @srrstats {ML3.0a} *As pre-processing, model specification, and training are controlled by the orsf() function, an input parameter (no_fit) enables models to be specified yet not fitted.*
#'
#' @param no_fit (_logical_) if `TRUE`, pre-processing steps are defined and
#'   parametrized, but training is not initiated. The object returned can be
#'   directly submitted to `orsf_train()` so long as `attach_data` is `TRUE`.
#'
#' @srrstats {ML3.0c} *when no_fit=TRUE, orsf() will return an object that can be directly trained using orsf_train().*
#'
#' and run the model training procedures.*
#' @param object an untrained aorsf object, created by setting
#'   `no_fit = TRUE` in `orsf()`.
#'
#' @return an accelerated oblique RSF object (`aorsf`)
#'
#' @details
#'
#' This function is based on and highly similar to the `ORSF` function
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
#' @srrstats {ML1.0} *Make a clear conceptual distinction between training and test data*
#'
#' __Training, out-of-bag error, and testing__
#'
#' In random forests, each tree is grown with a bootstrapped version of the training set. Because bootstrap samples are selected with replacement, each bootstrapped training set contains about two-thirds of instances in the original training set. The 'out-of-bag' data are instances that are _not_ in the bootstrapped training set. Each tree in the random forest can make predictions for its out-of-bag data, and the out-of-bag predictions can be aggregated to make an ensemble out-of-bag prediction. Since the out-of-bag data are not used to grow the tree, the accuracy of the ensemble out-of-bag predictions approximate the generalization error of the random forest. Generalization error refers to the error of a random forest's predictions when it is applied to predict outcomes for data that were not used to train it, i.e., testing data.
#'
#' __Missing data__
#' @srrstats {ML1.6a} *Explain why missing values are not admitted.*
#' Data passed to aorsf functions are not allowed to have missing values. A user should impute missing values using an R package with that purpose, such as `recipes` or `mlr3pipelines`. Other software such as `xgboost` send data with missing values down a decision tree based on whichever direction minimizes a specified error function. While this technique is very effective for axis-based decision trees, it is not clear how it should be applied in the case of oblique decision trees. For example, what should be done if three variables were used to split a node and one of these three variable has a missing value? In this case, mean imputation of the missing variable may be the best option.
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
#' @srrstats {G2.0a} *secondary documentation of arg lengths*
#' With the exception of `data_train` and `formula`, all inputs of `orsf()`
#'   should be an integer, double, or logical value of length 1.
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
#'
#' fit <- orsf(pbc_orsf, formula = Surv(time, status) ~ . - id)
#'
#' print(fit)
#'
#' @srrstats {ML1.6b} *Explicit example showing how missing values may be imputed rather than discarded.*
#'
#' \dontrun{requires too many external packages
#'
#' # --------------------------------------------------------------------------
#' # a standard machine learning workflow using aorsf and tidymodels
#' # --------------------------------------------------------------------------
#'
#' library(tidymodels)
#' library(tidyverse)
#' library(survivalROC)
#' library(aorsf)
#'
#' set.seed(329)
#'
#' # a recipe to impute missing values instead of discarding them.
#' # this is for illustration only. pbc_orsf does not have missing values
#'
#' imputer <- recipe(x = pbc_orsf, time + status ~ .) |>
#'  step_impute_mean(all_numeric_predictors()) |>
#'  step_impute_mode(all_nominal_predictors()) |>
#'  step_rm(id)
#'
#' # 10-fold cross validation; make a container for the pre-processed data
#' analyses <- vfold_cv(data = pbc_orsf, v = 10) |>
#'  mutate(recipe = map(splits, ~prep(imputer, training = training(.x))),
#'         data_train = map(recipe, juice),
#'         data_test = map2(splits, recipe, ~bake(.y, new_data = testing(.x))))
#'
#' # 10-fold cross validation; train models and compute test predictions
#' aorsf_data <- analyses |>
#'  select(data_train, data_test) |>
#'  mutate(fit = map(data_train, orsf, formula = time + status ~ .),
#'         pred = map2(fit, data_test, predict, pred_horizon = 3500),
#'         pred = map(pred, as.numeric))
#'
#' # testing sets are small, so pool them and compute 1 overall C-stat.
#' aorsf_eval <- aorsf_data |>
#'  select(data_test, pred) |>
#'  unnest(cols = everything()) |>
#'  summarize(
#'   auc = survivalROC(Stime = time,
#'                     status = status,
#'                     marker = pred,
#'                     predict.time = 3500,
#'                     span = 0.25*n()^(-0.20)) |>
#'    getElement('AUC')
#'  )
#'
#' # C-stat: 0.816
#' aorsf_eval
#'
#' # standard workflow for model development: fit and interpret
#'
#' aorsf_fit <- orsf(pbc_orsf, time + status ~ .,
#'                   importance = TRUE,
#'                   n_tree = 2500,
#'                   oobag_time = 3500)
#'
#' aorsf_smry <- orsf_summarize_uni(aorsf_fit, n_variables = 5)
#'
#' }
#'
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
                 oobag_fun = NULL,
                 importance = FALSE,
                 attach_data = TRUE,
                 no_fit = FALSE){

 #' @srrstats {G2.8} *As part of initial pre-processing, run checks on inputs to ensure that all other sub-functions receive inputs of a single defined class or type.*

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

 if(is.null(oobag_fun)){

  f_oobag_eval <- function(x) x
  type_oobag_eval <- 'H'

 } else {

  check_oobag_fun(oobag_fun)
  f_oobag_eval <- oobag_fun
  type_oobag_eval <- 'U'

 }

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


 types_y_data <- vector(mode = 'character', length = 2)

 for(i in seq_along(types_y_data)){
  types_y_data[i] <- class(data_train[[ names_y_data[i] ]])[1]
 }

 unit_y_names <- names_y_data[types_y_data == 'units']

 ui_y <- unit_info(data = data_train, .names = unit_y_names)

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

 #' @srrstats {G2.6} *ensure that one-dimensional inputs are appropriately pre-processed. aorsf does not deal with missing data as many other R packages are very good at dealing with it.*

 #' @srrstats {G2.13} *check for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*

 #' @srrstats {G2.15} *Never pass data with potential missing values to any base routines.*


 #' @srrstats {ML1.6} *do not admit missing values, and implement explicit pre-processing routines to identify whether data has any missing values. Throw errors appropriately and informatively when passed data contain missing values.*

 if(any(is.na(select_cols(data_train, c(names_y_data, names_x_data))))){

  stop("Please remove missing values from data_train, or impute them.",
       call. = FALSE)

 }

 #' @srrstats {G2.16} *Throw hard errors if undefined values are detected.*

 for(i in c(names_y_data, names_x_data)){

  if(any(is.infinite(data_train[[i]]))){
   stop("Please remove infinite values from ", i, ".",
        call. = FALSE)
  }

  if(any(is.nan(data_train[[i]]))){
   stop("Please remove NaN values from ", i, ".",
        call. = FALSE)
  }

 }

 fctr_check(data_train, names_x_data)
 fctr_id_check(data_train, names_x_data)

 fi <- fctr_info(data_train, names_x_data)

 #' @srrstats {G2.7} *aorsf accepts as input numeric and categorical predictor variables, including those with unit class. I do not think it is necessary to incorporate any other type of input, since it is relatively straightforward to convert data into a numeric or categorical format.*

 #' @srrstats {G2.11} *I cannot write code for every possible vector class to ensure that the vector data will be safely coerced into a a valid class and the attributes will be stored in the orsf_out object. It is much easier and safer for the user to convert a few columns to numeric or factor than it is for me to attempt writing code that will safely coerce every type of vector. That being said, I do find units columns to be helpful and I've written some code to make them an allowable class in input data.*

 types_x_data <- check_var_types(data_train,
                                 names_x_data,
                                 valid_types = c('numeric',
                                                 'integer',
                                                 'units',
                                                 'factor',
                                                 'ordered'))

 unit_x_names <- names_x_data[types_x_data == 'units']

 ui_x <- unit_info(data = data_train, .names = unit_x_names)

 names_x_numeric <- grep(pattern = "^integer$|^numeric$|^units$",
                         x = types_x_data)

 numeric_bounds <- NULL

 if(!is_empty(names_x_numeric)){
  numeric_bounds <-
   sapply(select_cols(data_train, names_x_data[names_x_numeric]),
          FUN = stats::quantile,
          probs = c(0.10, 0.25, 0.50, 0.75, 0.90))
 }

 y  <- as.matrix(select_cols(data_train, names_y_data))
 x  <- as.matrix(ref_code(data_train, fi, names_x_data))

 if(is.null(mtry)) mtry <- ceiling(sqrt(ncol(x)))

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

 #' @srrstats {ML2.3} *Values associated with transformations are recorded in the object returned by orsf(), specifically in object$forest[[<insert tree number>]]$x_mean*

 orsf_out <- orsf_fit(x                 = x_sort,
                      y                 = y_sort,
                      n_tree            = if(no_fit) 0 else n_tree,
                      n_split_          = n_split,
                      mtry_             = mtry,
                      leaf_min_events_  = leaf_min_events,
                      leaf_min_obs_     = leaf_min_obs,
                      split_min_events_ = split_min_events,
                      split_min_obs_    = split_min_obs,
                      cph_method_       = switch(tolower(cph_method),
                                                 'breslow' = 0,
                                                 'efron'   = 1),
                      cph_eps_          = cph_eps, #
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
                      type_beta_        = switch(orsf_type,
                                                 'cph' = 'C',
                                                 'net' = 'N'),
                      f_oobag_eval      = f_oobag_eval,
                      type_oobag_eval_  = type_oobag_eval)

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

  # clear labels for oobag evaluation type

  orsf_out$eval_oobag$stat_type <-
   switch(EXPR = orsf_out$eval_oobag$stat_type,
          'H' = "Harrell's C-statistic",
          'U' = "User-specified function")

  #' @srrstats {G2.10} *set drop = FALSE to ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behavior, and all column-extraction operations behave consistently regardless of the class of tabular data used as input.*

  orsf_out$surv_oobag <- orsf_out$surv_oobag[unsorted, , drop = FALSE]

 } else {

  # this would get added by orsf_fit if oobag_pred was TRUE
  orsf_out$pred_horizon <- stats::median(y[, 1])

 }

 n_leaves_mean <- 0

 if(!no_fit) {
  n_leaves_mean <-
   mean(sapply(orsf_out$forest, function(t) nrow(t$leaf_node_index)))
 }



 attr(orsf_out, 'mtry')               <- mtry
 attr(orsf_out, 'n_obs')              <- nrow(y_sort)
 attr(orsf_out, 'n_tree')             <- n_tree
 attr(orsf_out, 'names_y')            <- names_y_data
 attr(orsf_out, "names_x")            <- names_x_data
 attr(orsf_out, "names_x_ref")        <- colnames(x)
 attr(orsf_out, "types_x")            <- types_x_data
 attr(orsf_out, 'n_events')           <- n_events
 attr(orsf_out, 'max_time')           <- y_sort[nrow(y_sort), 1]
 attr(orsf_out, "unit_info")          <- c(ui_y, ui_x)
 attr(orsf_out, "fctr_info")          <- fi
 attr(orsf_out, 'n_leaves_mean')      <- n_leaves_mean
 attr(orsf_out, 'n_split')            <- n_split
 attr(orsf_out, 'leaf_min_events')    <- leaf_min_events
 attr(orsf_out, 'leaf_min_obs')       <- leaf_min_obs
 attr(orsf_out, 'split_min_events')   <- split_min_events
 attr(orsf_out, 'split_min_obs')      <- split_min_obs
 attr(orsf_out, 'cph_method')         <- cph_method
 attr(orsf_out, 'cph_eps')            <- cph_eps
 attr(orsf_out, 'cph_iter_max')       <- cph_iter_max
 attr(orsf_out, 'cph_pval_max')       <- cph_pval_max
 attr(orsf_out, 'cph_do_scale')       <- cph_do_scale
 attr(orsf_out, 'net_alpha')          <- net_alpha
 attr(orsf_out, 'net_df_target')      <- net_df_target
 attr(orsf_out, 'numeric_bounds')     <- numeric_bounds
 attr(orsf_out, 'trained')            <- !no_fit
 attr(orsf_out, 'n_retry')            <- n_retry
 attr(orsf_out, 'orsf_type')          <- orsf_type
 attr(orsf_out, 'f_beta')             <- f_beta
 attr(orsf_out, 'f_oobag_eval')       <- f_oobag_eval
 attr(orsf_out, 'type_oobag_eval')    <- type_oobag_eval
 attr(orsf_out, 'oobag_pred')         <- oobag_pred
 attr(orsf_out, 'oobag_eval_every')   <- oobag_eval_every
 attr(orsf_out, 'importance')         <- importance

 class(orsf_out) <- "aorsf"

 orsf_out


}

#' @srrstats {ML2.0a} *objects returned from orsf() with no_fit = TRUE can be directly submitted to orsf_train to train the model specification.*
#'
#' @rdname orsf
#' @export
orsf_train <- function(object){

 if(is_trained(object)){
  stop("object has already been trained", call. = FALSE)
 }

 if(is.null(object$data_train)){
  stop("object must have training data attached.",
       " Set attach_data = TRUE in orsf()",
       call. = FALSE)
 }

 y  <- as.matrix(select_cols(object$data_train, get_names_y(object)))

 x  <- as.matrix(ref_code(object$data_train,
                          get_fctr_info(object),
                          get_names_x(object, ref_code_names = FALSE)))

 sorted <- order(y[, 1],  # order this way for risk sets
                 -y[, 2]) # order this way for oob C-statistic.

 x_sort <- x[sorted, ]
 y_sort <- y[sorted, ]

 orsf_out <- orsf_fit(
  x                 = x_sort,
  y                 = y_sort,
  n_tree            = get_n_tree(object),
  n_split_          = get_n_split(object),
  mtry_             = get_mtry(object),
  leaf_min_events_  = get_leaf_min_events(object),
  leaf_min_obs_     = get_leaf_min_obs(object),
  split_min_events_ = get_split_min_events(object),
  split_min_obs_    = get_split_min_obs(object),
  cph_method_       = switch(tolower(get_cph_method(object)),
                             'breslow' = 0,
                             'efron'   = 1),
  cph_eps_          = get_cph_eps(object), #
  cph_iter_max_     = get_cph_iter_max(object),
  cph_pval_max_     = get_cph_pval_max(object),
  cph_do_scale_     = get_cph_do_scale(object),
  net_alpha_        = get_net_alpha(object),
  net_df_target_    = get_net_df_target(object),
  oobag_pred_       = get_oobag_pred(object),
  oobag_time_       = object$pred_horizon,
  oobag_eval_every_ = get_oobag_eval_every(object),
  oobag_importance_ = get_importance(object),
  max_retry_        = get_n_retry(object),
  f_beta            = get_f_beta(object),
  type_beta_        = switch(get_orsf_type(object),
                             'cph' = 'C',
                             'net' = 'N'),
  f_oobag_eval      = get_f_oobag_eval(object),
  type_oobag_eval_  = get_type_oobag_eval(object)
 )

 object$forest       <- orsf_out$forest
 object$surv_oobag   <- orsf_out$surv_oobag
 object$pred_horizon <- orsf_out$pred_horizon
 object$eval_oobag   <- orsf_out$eval_oobag
 object$importance   <- orsf_out$importance
 object$signif_means <- orsf_out$signif_means

 if(get_importance(object)){

  rownames(object$importance) <- colnames(x)

  object$importance <-
   rev(object$importance[order(object$importance), , drop=TRUE])

 }


 if(get_oobag_pred(object)){

  # put the oob predictions into the same order as the training data.
  unsorted <- vector(mode = 'integer', length = length(sorted))
  for(i in seq_along(unsorted)) unsorted[ sorted[i] ] <- i

  # clear labels for oobag evaluation type

  object$eval_oobag$stat_type <-
   switch(EXPR = object$eval_oobag$stat_type,
          'H' = "Harrell's C-statistic",
          'U' = "User-specified function")

  object$surv_oobag <- object$surv_oobag[unsorted, , drop = FALSE]

 }

 attr(object, "n_leaves_mean") <- mean(
  sapply(object$forest, function(t) nrow(t$leaf_node_index))
 )

 attr(object, 'trained') <- TRUE


 object

}









