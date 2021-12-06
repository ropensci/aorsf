

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
#' @param data (_data.frame_) that will be used to grow the forest.
#'
#' @param formula (_formula_) a formula object, with the response on the left
#'   of a `~` operator, and the terms on the right. See details.
#'
#' @param n_tree (_integer_) the number of trees to grow
#'
#' @param n_split (_integer_) the number of cut-points assessed when splitting
#'  a node in decision trees.
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
#' @param cph_method (_character_) a character string specifying the method
#'   for tie handling. If there are no tied death times all the methods are
#'   equivalent. Valid options are 'breslow' and 'efron'.
#'
#' @param cph_eps (_double_) When using Newton Raphson scoring to identify
#'   linear combinations of inputs, iteration continues in the algorithm
#'   until the relative change in  the log partial likelihood is less than
#'   `eps`, or the absolute change is less than `sqrt(eps)`. Must be positive.
#'
#' @param cph_iter_max (_integer_) When using Newton Raphson scoring to identify
#'   linear combinations of inputs, iteration continues until convergence
#'   (see `cph_eps` above) or the number of attempted iterations is equal to
#'   `cph_iter_max`.
#'
#' @param cph_pval_max (_double_) The maximum p-value allowed for a regression
#'   coefficient to remain non-zero. If the p-value for a given coefficient
#'   is above the maximum, the coefficient is set to zero and the variable
#'   no longer plays a role in the linear combination of inputs. Setting
#'   `cph_pval_max` to 1 ensures that every predict gets a non-zero
#'   coefficient in the linear combination of inputs.
#'
#' @param cph_do_scale (_logical_) if `TRUE`, values of predictors will be
#'   scaled prior to running Newton Raphson scoring. Setting to `FALSE` will
#'   reduce computation time but will also make the regression extremely
#'   unstable. Therefore, `orsf` will only let you set this input to `FALSE`
#'   if you also set `cph_iter_max` to 1.
#'
#' @param oobag_pred (_logical_) if `TRUE` out-of-bag predictions are returned
#'   in the `aorsf` object.
#'
#' @param oobag_eval_every (_integer_) The out-of-bag performance of the
#'   ensemble will be checked every `oobag_eval_every` trees. So, if
#'   `oobag_eval_every = 10`, then out-of-bag performance is checked
#'   after growing the 10th tree, the 20th tree, and so on.
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
#' fit <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 10)
#'
#'
orsf <- function(data_train,
                 formula,
                 n_tree = 500,
                 n_split = 5,
                 mtry = NULL,
                 leaf_min_events = 1,
                 leaf_min_obs = 15,
                 cph_method = 'breslow',
                 cph_eps = 1e-5,
                 cph_iter_max = 1,
                 cph_pval_max = 1,
                 cph_do_scale = TRUE,
                 oobag_pred = FALSE,
                 oobag_time = NULL,
                 oobag_eval_every = n_tree,
                 importance = FALSE,
                 attach_data = TRUE){

 # Run checks
 Call <- match.call()

 check_call(
  Call,
  expected = list(
   'data_train' = list(
    class = 'data.frame'
   ),
   'formula' = list(
    class = 'formula'
   ),
   'n_tree' = list(
    type = 'numeric',
    length = 1,
    lwr = 1,
    integer = TRUE
   ),
   'n_split' = list(
    type = 'numeric',
    length = 1,
    lwr = 1,
    integer = TRUE
   ),
   'leaf_min_events' = list(
    type = 'numeric',
    length = 1,
    lwr = 1,
    integer = TRUE
   ),
   'leaf_min_obs' = list(
    type = 'numeric',
    length = 1,
    lwr = 1,
    integer = TRUE
   ),
   'cph_method' = list(
    type = 'character',
    options = c("breslow", "efron")
   ),
   'cph_eps' = list(
    type = 'numeric',
    length = 1,
    lwr = 0,
    integer = FALSE
   ),
   'cph_iter_max' = list(
    type = 'numeric',
    length = 1,
    lwr = 1,
    integer = TRUE
   ),
   'cph_pval_max' = list(
    type = 'numeric',
    length = 1,
    lwr = 0,
    upr = 1
   ),
   'cph_do_scale' = list(
    type = 'logical',
    length = 1
   ),
   'oobag_pred' = list(
    type = 'logical',
    length = 1
   ),
   oobag_time = list(
    type = 'numeric',
    length = 1,
    lwr = 0
   ),
   'oobag_eval_every' = list(
    type = 'numeric',
    integer = TRUE,
    lwr = 1,
    upr = n_tree
   ),
   'attach_data' = list(
    type = 'logical',
    length = 1
   )
  )
 )

 if(importance && !oobag_pred) oobag_pred <- TRUE # Should I add a warning?

 formula_terms <- stats::terms(formula, data=data_train)

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
   "variables in formula were not found in ", deparse(Call$data_train),
   " data_train: ",
   paste_collapse(names_not_found, last = ' and ')
  )
  stop(msg, call. = FALSE)
 }

 if(any(is.na(data_train[, c(names_y_data, names_x_data)]))){
  stop("Please remove missing values from ",
       deparse(Call$data_train),
       " or impute them",
       call. = FALSE)
 }

 fctr_check(data_train, names_x_data)
 fi <- fctr_info(data_train, names_x_data)
 y  <- as.matrix(data_train[, names_y_data])
 x  <- as.matrix(one_hot(data_train, fi, names_x_data))

 types_x_data <- check_var_types(data_train,
                                 names_x_data,
                                 valid_types = c('numeric',
                                                 'integer',
                                                 'factor',
                                                 'ordered'))

 check_arg_uni(arg_value = y[,2],
               arg_name = paste0(deparse(Call$data_train),
                                 '$', names_y_data[2]),
               expected_uni = c(0,1))

 check_arg_is(arg_value = y[,1],
              arg_name = paste0(deparse(Call$data_train),
                                '$', names_y_data[1]),
              expected_class = 'numeric')

 check_arg_gt(arg_value = y[,1],
              arg_name = paste0(deparse(Call$data_train),
                                '$', names_y_data[1]),
              bound = 0)

 if(is.null(mtry)){

  mtry <- ceiling(sqrt(ncol(x)))

 } else {

  check_arg_is_integer(arg_name = 'mtry',
                       arg_value = mtry)

  check_bound_lwr(arg_name = 'mtry',
                  arg_value = mtry,
                  bound_lwr = 1)

  check_arg_length(arg_name = 'mtry',
                   arg_value = mtry,
                   expected_length = 1)

 }

 if(!is.null(oobag_time)){
  if(oobag_time == 0)

   stop("Out of bag prediction time (oobag_time) must be > 0",
        call. = FALSE)

 } else {

  # sneaky way to tell orsf.cpp to make up its own oobag_time
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
                      oobag_pred_       = oobag_pred,
                      oobag_time_       = oobag_time,
                      oobag_eval_every_ = oobag_eval_every,
                      oobag_importance_ = importance)

 orsf_out$data_train <- if(attach_data) data_train else NULL

 if(importance){
  rownames(orsf_out$importance) <- colnames(x)
 }

 names_x_numeric <- grep(pattern = "^integer$|^numeric$",
                         x = types_x_data)

 numeric_bounds <- NULL

 if(!is_empty(names_x_numeric)){
  numeric_bounds <- sapply(data_train[, names_x_data[names_x_numeric] ],
                           quantile, probs = c(0.10, 0.20, 0.80, 0.90))
 }

 if(oobag_pred){

  # put the oob predictions into the same order as the training data.
  unsorted <- vector(mode = 'integer', length = length(sorted))
  for(i in seq_along(unsorted)) unsorted[ sorted[i] ] <- i

  orsf_out$surv_oobag <- orsf_out$surv_oobag[unsorted, , drop = FALSE]

 }

 n_leaves_mean <-
  mean(sapply(orsf_out$forest, function(t) nrow(t$leaf_node_index)))

 class(orsf_out) <- "aorsf"

 attr(orsf_out, 'mtry')            <- mtry
 attr(orsf_out, 'n_obs')           <- nrow(y_sort)
 attr(orsf_out, 'n_tree')          <- n_tree
 attr(orsf_out, "names_x")         <- names_x_data
 attr(orsf_out, "types_x")         <- types_x_data
 attr(orsf_out, 'n_events')        <- sum(y_sort[,2])
 attr(orsf_out, 'max_time')        <- y_sort[nrow(y_sort), 1]
 attr(orsf_out, "fctr_info")       <- fi
 attr(orsf_out, 'n_leaves_mean')   <- n_leaves_mean
 attr(orsf_out, 'n_split')         <- n_split
 attr(orsf_out, 'leaf_min_events') <- leaf_min_events
 attr(orsf_out, 'leaf_min_obs')    <- leaf_min_obs
 attr(orsf_out, 'numeric_bounds')  <- numeric_bounds

 orsf_out


}
