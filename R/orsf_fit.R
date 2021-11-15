


#' Fit oblique RSF
#'
#' The oblique random survival forest (RSF) is an extension of the RSF
#'   algorithm developed by Ishwaran et al and maintained in the
#'   `RandomForestSRC` package. The only difference between oblique
#'   RSF and Ishwaran's RSF is that oblique RSFs use linear combinations
#'   of input variables instead of using the input variable as-is when
#'   growing new nodes in survival decision trees. For more details on
#'   the oblique RSF, see Jaeger et al, 2019.
#'
#' This function is based on and highly similar to the `ORSF` function
#'   in the `obliqueRSF` R package. The primary difference is that this
#'   function runs about 100 times faster because it uses a simplified
#'   Newton Raphson scoring algorithm to identify linear combinations of
#'   inputs rather than performing penalized regression using routines in
#'   `glmnet`.
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
#' @param cph_pval_max (_double_)
#'
#' @param oobag_pred (_logical_) if `TRUE` out-of-bag predictions are returned
#'   in the `aorsf` object.
#'
#' @return an accelerated oblique RSF object (`aorsf`)
#'
#' @details
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
orsf <- function(data,
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
                 oobag_pred = FALSE){

 # Run checks
 Call <- match.call()

 check_call(
  Call,
  expected = list(
   'data' = list(
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
    lwr = 1,
    integer = TRUE
   ),
   'oobag_pred' = list(
    type = 'logical',
    length = 1
   )
  )
 )

 formula_terms <- stats::terms(formula, data=data)

 if(attr(formula_terms, 'response') == 0)
  stop("formula must have a response", call. = FALSE)

 if(length(attr(formula_terms, 'term.labels')) < 2)
  stop("formula must have at least 2 predictors", call. = FALSE)

 names_y_data <- all.vars(formula[[2]])

 if(length(names_y_data) != 2)
  stop("formula must have two variables (time & status) as the response",
       call. = FALSE)

 names_x_data <- attr(formula_terms, 'term.labels')

 names_not_found <- setdiff(c(names_y_data, names_x_data), names(data))

 if(!is_empty(names_not_found)){
  msg <- paste0(
   "variables in formula were not found in ", deparse(Call$data),
   " data: ",
   paste_collapse(names_not_found, last = ' and ')
  )
  stop(msg, call. = FALSE)
 }

 if(any(is.na(data[, c(names_y_data, names_x_data)]))){
  stop("Please remove missing values from ",
       deparse(Call$data),
       " or impute them",
       call. = FALSE)
 }

 fctr_check(data, names_x_data)
 fi <- fctr_info(data, names_x_data)
 y  <- as.matrix(data[, names_y_data])
 x  <- as.matrix(one_hot(data, fi, names_x_data))

 types_x_data <- check_var_types(data,
                                 names_x_data,
                                 valid_types = c('numeric',
                                                 'integer',
                                                 'factor',
                                                 'ordered'))

 check_arg_uni(arg_value = y[,2],
               arg_name = paste0(deparse(Call$data), '$', names_y_data[2]),
               expected_uni = c(0,1))

 check_arg_is(arg_value = y[,1],
              arg_name = paste0(deparse(Call$data), '$', names_y_data[1]),
              expected_class = 'numeric')

 check_arg_gt(arg_value = y[,1],
              arg_name = paste0(deparse(Call$data), '$', names_y_data[1]),
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

 sorted <- order(y[, 1])

 x_sort <- x[sorted, ]
 y_sort <- y[sorted, ]


 orsf_out <- orsf_fit(x                = x_sort,
                      y                = y_sort,
                      n_tree           = n_tree,
                      n_split_         = n_split,
                      mtry_            = mtry,
                      leaf_min_events_ = leaf_min_events,
                      leaf_min_obs_    = leaf_min_obs,
                      cph_method_      = switch(tolower(cph_method),
                                                'breslow' = 0,
                                                'efron'   = 1),
                      cph_eps_         = cph_eps,
                      cph_iter_max_    = cph_iter_max,
                      cph_pval_max_    = cph_pval_max,
                      oobag_pred_      = oobag_pred)

 orsf_out$fctr_info <- fi
 orsf_out$names_x <- names_x_data
 orsf_out$types_x <- types_x_data

 class(orsf_out) <- "aorsf"

 attr(orsf_out, 'max_time') = y_sort[nrow(y_sort), 1]

 return(orsf_out)


}
