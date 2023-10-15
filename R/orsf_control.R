


#' Accelerated ORSF control
#'
#' Fast methods to identify linear combinations of predictors while
#'   fitting an [orsf] model.
#'
#' @param method (_character_) a character string specifying the method
#'   for tie handling. If there are no ties, all the methods are
#'   equivalent. Valid options are 'breslow' and 'efron'. The Efron
#'   approximation is the default because it is more accurate when dealing
#'   with tied event times and has similar computational efficiency compared
#'   to the Breslow method.
#'
#' @param do_scale (_logical_) if `TRUE`, values of predictors will be
#'   scaled prior to each instance of Newton Raphson scoring, using summary
#'   values from the data in the current node of the decision tree.
#'
#' @param ... `r roxy_dots()`
#'
#' @return an object of class `'orsf_control'`, which should be used as
#'  an input for the `control` argument of [orsf].
#'
#' @export
#'
#' @family orsf_control
#'
#' @details
#'
#'  code from the  [survival package](https://github.com/therneau/survival/blob/master/src/coxfit6.c)
#'   was modified to make this routine.
#'
#' Adjust `do_scale` _at your own risk_. Setting `do_scale = FALSE` will
#'  reduce computation time but will also make the `orsf` model dependent
#'  on the scale of your data, which is why the default value is `TRUE`.
#'
#'
#' @examples
#'
#' orsf(data = pbc_orsf,
#'      formula = Surv(time, status) ~ . - id,
#'      control = orsf_control_fast())
#'
orsf_control_fast <- function(method = 'efron',
                              do_scale = TRUE,
                              ...){

 check_dots(list(...), orsf_control_fast)

 method <- tolower(method)

 check_control_cph(method = method, do_scale = do_scale)

 ties_method <- method

 orsf_control(tree_type = 'unknown',
              method = 'glm',
              scale_x = do_scale,
              ties = ties_method,
              net_mix = 0.5,
              target_df = NULL,
              max_iter = 1,
              epsilon = 1e-9)

}


#' Cox regression ORSF control
#'
#' Use the coefficients from a proportional hazards model
#'  to create linear combinations of predictor variables
#'  while fitting an [orsf] model.
#'
#' @inheritParams orsf_control_fast
#'
#' @param eps (_double_) When using Newton Raphson scoring to identify
#'   linear combinations of inputs, iteration continues in the algorithm
#'   until the relative change in  the log partial likelihood is less than
#'   `eps`, or the absolute change is less than `sqrt(eps)`. Must be positive.
#'   A default value of 1e-09 is used for consistency with
#'   [survival::coxph.control].
#'
#' @param iter_max (_integer_) iteration continues until convergence
#'   (see `eps` above) or the number of attempted iterations is equal to
#'   `iter_max`.
#'
#' @return an object of class `'orsf_control'`, which should be used as
#'  an input for the `control` argument of [orsf].
#'
#' @export
#'
#' @family orsf_control
#'
#' @details
#'
#'  code from the  [survival package](https://github.com/therneau/survival/blob/master/src/coxfit6.c)
#'   was modified to make this routine.
#'
#'  For more details on the Cox proportional hazards model, see
#'  [coxph][survival::coxph] and/or Therneau and Grambsch (2000).
#'
#' @references
#'
#' Therneau T.M., Grambsch P.M. (2000) The Cox Model. In: Modeling Survival
#'   Data: Extending the Cox Model. Statistics for Biology and Health.
#'   Springer, New York, NY. DOI: 10.1007/978-1-4757-3294-8_3
#'
#' @examples
#'
#' orsf(data = pbc_orsf,
#'      formula = Surv(time, status) ~ . - id,
#'      control = orsf_control_cph())
#'
orsf_control_cph <- function(method = 'efron',
                             eps = 1e-9,
                             iter_max = 20,
                             ...){


 method <- tolower(method)

 check_dots(list(...), orsf_control_cph)

 check_control_cph(method = method,
                   eps = eps,
                   iter_max = iter_max)

 ties_method <- method

 orsf_control(tree_type = 'unknown',
              method = 'glm',
              scale_x = TRUE,
              ties = ties_method,
              net_mix = 0.5,
              target_df = NULL,
              max_iter = iter_max,
              epsilon = eps)

}

#' Penalized Cox regression ORSF control
#'
#' Use regularized Cox proportional hazard models to identify linear
#'   combinations of input variables while fitting an [orsf] model.
#'
#' @param alpha (_double_) The elastic net mixing parameter. A value of 1 gives the
#'  lasso penalty, and a value of 0 gives the ridge penalty. If multiple
#'  values of alpha are given, then a penalized model is fit using each
#'  alpha value prior to splitting a node.
#'
#' @param df_target (_integer_) Preferred number of variables used in a linear combination.
#'
#' @param ... `r roxy_dots()`
#'
#' @inherit orsf_control_cph return
#'
#' @details
#'
#' `df_target` has to be less than `mtry`, which is a separate argument in
#'  [orsf] that indicates the number of variables chosen at random prior to
#'  finding a linear combination of those variables.
#'
#' @export
#'
#' @family orsf_control
#'
#' @references
#'
#' `r roxy_cite_simon_2011()`
#'
#'
#'
#' @examples
#'
#' # orsf_control_net() is considerably slower than orsf_control_cph(),
#' # The example uses n_tree = 25 so that my examples run faster,
#' # but you should use at least 500 trees in applied settings.
#'
#' orsf(data = pbc_orsf,
#'      formula = Surv(time, status) ~ . - id,
#'      n_tree = 25,
#'      control = orsf_control_net())

orsf_control_net <- function(alpha = 1/2,
                             df_target = NULL,
                             ...){

 check_dots(list(...), orsf_control_net)
 check_control_net(alpha, df_target)

 orsf_control(tree_type = 'unknown',
              method = 'net',
              scale_x = TRUE,
              ties = 'efron',
              net_mix = alpha,
              target_df = df_target,
              max_iter = 20,
              epsilon = 1e-9)

}

#' Custom ORSF control
#'
#' @param beta_fun (_function_) a function to define coefficients used
#'  in linear combinations of predictor variables. `beta_fun` must accept
#'  three inputs named `x_node`, `y_node` and `w_node`, and should expect
#'  the following types and dimensions:
#'
#' - `x_node` (_matrix_; `n` rows, `p` columns)
#' - `y_node` (_matrix_; `n` rows, `2` columns)
#' - `w_node` (_matrix_; `n` rows, `1` column)
#'
#' In addition, `beta_fun` must return a matrix with p rows and 1 column. If
#'  any of these conditions are not met, `orsf_control_custom()` will let
#'  you know.
#'
#' @param ... `r roxy_dots()`
#'
#' @inherit orsf_control_cph return
#'
#' @export
#'
#' @family orsf_control
#'
#' @includeRmd Rmd/orsf_control_custom_examples.Rmd

orsf_control_custom <- function(beta_fun, ...){

 check_dots(list(...), .f = orsf_control_custom)
 check_beta_fun(beta_fun)

 orsf_control(tree_type = 'unknown',
              method = beta_fun,
              scale_x = TRUE,
              ties = 'efron',
              net_mix = 0.5,
              target_df = NULL,
              max_iter = 20,
              epsilon = 1e-9)


}


#' Oblique random forest control
#'
#' @param tree_type (_character_) the type of tree. Valid options are
#'
#'  - "classification", i.e., categorical outcomes
#'  - "regression", i.e., continuous outcomes
#'  - "survival", i.e., time-to event outcomes
#'
#' @param method (_character_ or _function_) how to identify linear
#'  linear combinations of predictors. If `method` is a character value,
#'  it must be one of:
#'
#'  - 'glm': linear, logistic, and cox regression
#'  - 'net': same as 'glm' but with penalty terms
#'  - 'pca': principal component analysis
#'  - 'random': random draw from uniform distribution
#'
#' If `method` is a _function_, it will be used to identify  linear
#'  combinations of predictor variables. `method` must in this case accept
#'  three inputs named `x_node`, `y_node` and `w_node`, and should expect
#'  the following types and dimensions:
#'
#'  - `x_node` (_matrix_; `n` rows, `p` columns)
#'  - `y_node` (_matrix_; `n` rows, `2` columns)
#'  - `w_node` (_matrix_; `n` rows, `1` column)
#'
#' In addition, `method` must return a matrix with p rows and 1 column.
#'
#' @param scale_x (_logical_) if `TRUE`, values of predictors will be
#'   scaled prior to each instance of finding a linear combination of
#'   predictors, using summary values from the data in the current node
#'   of the decision tree.
#'
#' @param ties (_character_) a character string specifying the method
#'   for tie handling. Only relevant when modeling survival outcomes
#'   and using a method that engages with tied outcome times.
#'   If there are no ties, all the methods are equivalent. Valid options
#'   are 'breslow' and 'efron'. The Efron approximation is the default
#'   because it is more accurate when dealing with tied event times and
#'   has similar computational efficiency compared to the Breslow method.
#'
#' @param net_mix (_double_) The elastic net mixing parameter. A value of 1
#'  gives the lasso penalty, and a value of 0 gives the ridge penalty. If
#'  multiple values of alpha are given, then a penalized model is fit using
#'  each alpha value prior to splitting a node.
#'
#' @param target_df (_integer_) Preferred number of variables used in each
#'   linear combination. For example, with `mtry` of 5 and `target_df` of 3,
#'   we sample 5 predictors and look for the best linear combination using
#'   3 of them.
#'
#' @param max_iter  (_integer_) iteration continues until convergence
#'   (see `eps` above) or the number of attempted iterations is equal to
#'   `iter_max`.
#'
#' @param epsilon (_double_) When using most modeling based method,
#'   iteration continues in the algorithm until the relative change in
#'   some kind of objective is less than `epsilon`, or the absolute
#'   change is less than `sqrt(epsilon)`.
#'
#' @param ... `r roxy_dots()`
#'
#' @details
#'
#' Adjust `scale_x` _at your own risk_. Setting `scale_x = FALSE` will
#'  reduce computation time but will also make the `orsf` model dependent
#'  on the scale of your data, which is why the default value is `TRUE`.
#'
#'
#' @return an object of class `'orsf_control'`, which should be used as
#'  an input for the `control` argument of [orsf]. Components are:
#'
#' - `tree_type`: type of trees to fit
#' - `lincomb_type`: method for linear combinations
#' - `lincomb_eps`: epsilon for convergence
#' - `lincomb_iter_max`: max iterations
#' - `lincomb_scale`: to scale or not.
#' - `lincomb_alpha`: mixing parameter
#' - `lincomb_df_target`: target degrees of freedom
#' - `lincomb_ties_method`: method for ties in survival time
#' - `lincomb_R_function`: R function for custom splits
#'
#' @examples
#'
#' orsf_control_classification()
#' orsf_control_regression()
#' orsf_control_survival()
#'
orsf_control <- function(tree_type,
                         method,
                         scale_x,
                         ties,
                         net_mix,
                         target_df,
                         max_iter,
                         epsilon,
                         ...){

 custom <- is.function(method)

 if(custom){
  lincomb_R_function <- method
 } else if (method == 'net') {
  lincomb_R_function <- penalized_cph
 } else {
  lincomb_R_function <- function(x) x
 }

 structure(
  .Data = list(
   tree_type = tree_type,
   lincomb_type = if(custom) "custom" else method,
   lincomb_eps = epsilon,
   lincomb_iter_max = max_iter,
   lincomb_scale = scale_x,
   lincomb_alpha = net_mix,
   lincomb_df_target = target_df,
   lincomb_ties_method = ties,
   lincomb_R_function = lincomb_R_function
  ),
  class = c(paste('orsf_control', tree_type, sep = '_'),
            'orsf_control')
 )

}

#' @rdname orsf_control
#' @export
orsf_control_classification <- function(method = 'glm',
                                        scale_x = TRUE,
                                        net_mix = 0.5,
                                        target_df = NULL,
                                        max_iter = 20,
                                        epsilon = 1e-9,
                                        ...){

 check_dots(list(...), orsf_control_classification)

 orsf_control("classification",
              method = method,
              scale_x = scale_x,
              ties = 'efron',
              net_mix = net_mix,
              target_df = target_df,
              max_iter = max_iter,
              epsilon = epsilon,
              ...)

}

#' @rdname orsf_control
#' @export
orsf_control_regression <- function(method = 'glm',
                                    scale_x = TRUE,
                                    net_mix = 0.5,
                                    target_df = NULL,
                                    max_iter = 20,
                                    epsilon = 1e-9,
                                    ...){

 check_dots(list(...), orsf_control_regression)

 orsf_control("regression",
              method = method,
              scale_x = scale_x,
              ties = 'efron',
              net_mix = net_mix,
              target_df = target_df,
              max_iter = max_iter,
              epsilon = epsilon,
              ...)

}

#' @rdname orsf_control
#' @export
orsf_control_survival <- function(method = 'glm',
                                  scale_x = TRUE,
                                  ties = 'efron',
                                  net_mix = 0.5,
                                  target_df = NULL,
                                  max_iter = 20,
                                  epsilon = 1e-9,
                                  ...){

 check_dots(list(...), orsf_control_survival)

 orsf_control("survival",
              method = method,
              scale_x = scale_x,
              ties = ties,
              net_mix = net_mix,
              target_df = target_df,
              max_iter = max_iter,
              epsilon = epsilon,
              ...)

}

