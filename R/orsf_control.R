


#' Accelerated ORSF control
#'
#' @srrstats {ML3.5} *The orsf_control_ function family allows users to control optimization algorithms used to grow random forests*
#' @srrstats {G1.4} *documented with Roxygen*
#' @srrstats {ML2.4} *Default values of all transformations are explicitly documented.*
#' @srrstats {G1.3} *clarify Newton-Raphson scoring and Cox PH.*
#' @srrstats {ML3.5a} *Specify Newton-Raphson scoring or penalized regression as the type of algorithm used to explore the search space, i.e., the space of possible linear combinations*
#' @srrstats {ML3.6, ML3.6a} *Implement usage of multiple ways of exploring search space*
#' @srrstats {G3.0} *use eps to avoid comparing floating point numbers for equality*
#' @srrstats {ML2.5} *Provide the option to bypass default transformations.*
#'
#' Use a single iteration of Newton-Raphson scoring to identify linear
#' combinations of predictors while fitting an [orsf] model.
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
#' @param ... Further arguments passed to or from other methods
#'   (not currently used).
#'
#' @return an object of class `'aorsf_control'`, which should be used as
#'  an input for the `control` argument of [orsf].
#'
#' @export
#'
#' @family orsf_control
#'
#' @details
#'
#'  For more details on of Newton-Raphson scoring and the Cox proportional
#'  hazards model, see Therneau and Grambsch (2000).
#'
#' Adjust `do_scale` _at your own risk_. Setting `do_scale = FALSE` will
#'  reduce computation time but will also make the `orsf` model dependent
#'  on the scale of your data, which is why the default value is `TRUE`. It
#'  would be a good idea to center and scale your predictors prior to running
#'  `orsf()` if you plan on setting `do_scale = FALSE`.
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
#'      control = orsf_control_fast())
#'
orsf_control_fast <- function(method = 'efron',
                              do_scale = TRUE,
                              ...){

 check_dots(list(...), orsf_control_fast)

 check_control_cph(method = method, do_scale = do_scale)

 structure(
  .Data = list(cph_method = method,
               cph_eps = 1e-9,
               cph_iter_max = 1,
               cph_do_scale = do_scale),
  class = 'aorsf_control',
  type = 'fast'
 )

}


#' Cox proportional hazards ORSF control
#'
#'
#' Use the coefficients from a proportional hazards model
#'  to create linear combinations of predictor variables
#'  while fitting an [orsf] model.
#'
#' @param method (_character_) a character string specifying the method
#'   for tie handling. If there are no ties, all the methods are
#'   equivalent. Valid options are 'breslow' and 'efron'. The Efron
#'   approximation is the default because it is more accurate when dealing
#'   with tied event times and has similar computational efficiency compared
#'   to the Breslow method.
#'
#' @param eps (_double_) When using Newton Raphson scoring to identify
#'   linear combinations of inputs, iteration continues in the algorithm
#'   until the relative change in  the log partial likelihood is less than
#'   `eps`, or the absolute change is less than `sqrt(eps)`. Must be positive.
#'   A default value of 1e-09 is used for consistency with
#'   [survival::coxph.control].
#'
#' @param iter_max (_integer_) When using Newton Raphson scoring to identify
#'   linear combinations of inputs, iteration continues until convergence
#'   (see `eps` above) or the number of attempted iterations is equal to
#'   `iter_max`.
#'
#' @param ... Further arguments passed to or from other methods
#'   (not currently used).
#'
#' @return an object of class `'aorsf_control'`, which should be used as
#'  an input for the `control` argument of [orsf].
#'
#' @export
#'
#' @family orsf_control
#'
#' @details
#'
#'  For more details on of Newton-Raphson scoring and the Cox proportional
#'  hazards model, see Therneau and Grambsch (2000).
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


 #' @srrstats {G2.3b} *ensure input of character parameters is not case dependent*
 method <- tolower(method)

 check_dots(list(...), orsf_control_cph)

 check_control_cph(method = method,
                   eps = eps,
                   iter_max = iter_max)

 structure(
  .Data = list(cph_method = method,
               cph_eps = eps,
               cph_iter_max = iter_max,
               cph_do_scale = TRUE),
  class = 'aorsf_control',
  type = 'cph'
 )

}

#' Elastic net control
#'
#' @srrstats {ML3.5a} *Specify regularization of the coxph model as the type of algorithm used to explore the search space*
#'
#' @srrstats {ML3.6b} *Use the loss function associated with the penalized coxph model instead of the Newton Raphson scoring algorithm (default).*

#' Use regularized Cox proportional hazard models to identify linear
#'   combinations of input variables while fitting an [orsf] model.
#'
#' @param alpha The elastic net mixing parameter. A value of 1 gives the
#'  lasso penalty, and a value of 0 gives the ridge penalty. If multiple
#'  values of alpha are given, then a penalized model is fit using each
#'  alpha value prior to splitting a node.
#'
#' @param df_target Preferred number of variables used in a linear combination.
#'  Note: this has to be less than `mtry`, which is a separate argument in
#'  [orsf] that indicates the number of variables chosen at random prior to
#'  finding a linear combination of those variables.
#'
#' @param ... Further arguments passed to or from other methods
#'   (not currently used).
#'
#' @inherit orsf_control_cph return
#'
#' @export
#'
#' @family orsf_control
#'
#' @references
#'
#'  Simon N, Friedman J, Hastie T, Tibshirani R. Regularization paths for Cox’s proportional hazards model via coordinate descent. *Journal of statistical software*. 2011 Mar;39(5):1. DOI: 10.18637/jss.v039.i05
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

 structure(
  .Data = list(net_alpha = alpha,
               net_df_target = df_target),
  class = 'aorsf_control',
  type = 'net'
 )

}

#' Custom control of oblique decision trees
#'
#' @param beta_fun (_function_) a function to define coefficients used
#'  in linear combinations of predictor variables. `beta_fun` must accept
#'  three inputs named `x_node`, `y_node` and `w_node`, and should expect
#'  the following types and dimensions:
#'
#' - `x_node` (_matrix_; n rows, p columns)
#' - `y_node` (_matrix_; n rows, 2 columns)
#' - `w_node` (_matrix_; n rows, 1 column)
#'
#' In addition, `beta_fun` must return a matrix with p rows and 1 column. If
#'  any of these conditions are not met, `orsf_control_custom()` will let
#'  you know.
#'
#' @param ... Further arguments passed to or from other methods
#'   (not currently used).
#'
#' @inherit orsf_control_cph return
#'
#' @export
#'
#' @family orsf_control
#'
#' @examples
#'
#' # fit an oblique random survival forest using random coefficients to
#' # generate linear combinations of predictor variables. First, define
#' # a function that supplies the random coefficients:
#'
#' f <- function(x_node, y_node, w_node) { matrix(runif(ncol(x_node)), ncol=1) }
#'
#' # next, plug the function into orsf_control_custom(), which is in turn
#' # passed into orsf():
#'
#' fit_rando <- orsf(pbc_orsf,
#'                   Surv(time, status) ~ .,
#'                   control = orsf_control_custom(beta_fun = f),
#'                   n_tree = 500)
#'
#' # last, check the out-of-bag performance.
#' # it's surprising how well the random approach works.
#'
#' fit_rando$eval_oobag

orsf_control_custom <- function(beta_fun, ...){

 check_dots(list(...), .f = orsf_control_custom)
 check_beta_fun(beta_fun)

 structure(
  .Data = list(beta_fun = beta_fun),
  class = 'aorsf_control',
  type = 'custom'
 )


}



