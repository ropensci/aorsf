

#' Cox proportional hazards control
#'
#' @srrstats {G1.4} *documented with Roxygen*
#'
#' @srrstats {G1.3} clarify Newton-Raphson scoring and Cox PH.
#'
#' Use Newton-Raphson scoring to identify linear combinations of input
#'   variables while fitting an [orsf] model. For more details on
#'   of Newton-Raphson scoring and the Cox proportional hazards
#'   model, see Therneau and Grambsch (2000).
#'
#' @param method (_character_) a character string specifying the method
#'   for tie handling. If there are no ties, all the methods are
#'   equivalent. Valid options are 'breslow' and 'efron'.
#'
#' @param eps (_double_) When using Newton Raphson scoring to identify
#'   linear combinations of inputs, iteration continues in the algorithm
#'   until the relative change in  the log partial likelihood is less than
#'   `eps`, or the absolute change is less than `sqrt(eps)`. Must be positive.
#'
#' @param iter_max (_integer_) When using Newton Raphson scoring to identify
#'   linear combinations of inputs, iteration continues until convergence
#'   (see `eps` above) or the number of attempted iterations is equal to
#'   `iter_max`.
#'
#' @param pval_max (_double_) The maximum p-value allowed for a regression
#'   coefficient to remain non-zero. If the p-value for a given coefficient
#'   is above the maximum, the coefficient is set to zero and the variable
#'   no longer plays a role in the linear combination of inputs. Setting
#'   `pval_max` to 1 ensures that every predict gets a non-zero
#'   coefficient in the linear combination of inputs.
#'
#' @param do_scale (_logical_) if `TRUE`, values of predictors will be
#'   scaled prior to running Newton Raphson scoring. Setting to `FALSE` will
#'   reduce computation time but will also make the regression extremely
#'   unstable. Therefore, `orsf` will only let you set this input to `FALSE`
#'   if you also set `iter_max` to 1.
#'
#' @return an object of class `'aorsf_control'`, which should be used as
#'  an input for the `control` argument of [orsf].
#'
#' @export
#'
#' @references
#'
#' Therneau T.M., Grambsch P.M. (2000) The Cox Model. In: Modeling Survival
#'   Data: Extending the Cox Model. Statistics for Biology and Health.
#'   Springer, New York, NY. DOI: 10.1007/978-1-4757-3294-8_3
#'
#' @examples
#'
#' orsf(data_train = pbc_orsf,
#'      formula = Surv(time, status) ~ . - id,
#'      control = orsf_control_cph())
#'
orsf_control_cph <- function(method = 'breslow',
                             eps = 1e-5,
                             iter_max = 1,
                             pval_max = 1,
                             do_scale = TRUE){


 check_control_cph(method, eps, iter_max, pval_max, do_scale)

 structure(
  .Data = list(cph_method = method,
               cph_eps = eps,
               cph_iter_max = iter_max,
               cph_pval_max = pval_max,
               cph_do_scale = do_scale),
  class = 'aorsf_control',
  type = 'cph'
 )

}

#' Elastic net control
#'
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
#' @return an object of class `'aorsf_control'`, which should be used as
#'  an input for the `control` argument of [orsf].
#'
#' @export
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
#' orsf(data_train = pbc_orsf,
#'      formula = Surv(time, status) ~ . - id,
#'      n_tree = 25,
#'      control = orsf_control_net())

orsf_control_net <- function(alpha = 1/2,
                             df_target = NULL){

 check_control_net(alpha, df_target)

 structure(
  .Data = list(net_alpha = alpha,
               net_df_target = df_target),
  class = 'aorsf_control',
  type = 'net'
 )

}



