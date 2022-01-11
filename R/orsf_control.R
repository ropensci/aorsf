

#' Cox proportional hazards control
#'
#' @param method (_character_) a character string specifying the method
#'   for tie handling. If there are no tied death times all the methods are
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

check_control_cph <- function(method, eps, iter_max, pval_max, do_scale){

 check_arg_type(arg_value = method,
                arg_name = 'method',
                expected_type = 'character')

 check_arg_is_valid(arg_value = method,
                    arg_name = 'method',
                    valid_options = c("breslow", "efron"))


 check_arg_type(arg_value = eps,
                arg_name = 'eps',
                expected_type = 'numeric')

 check_arg_gt(arg_value = eps,
              arg_name = 'eps',
              bound = 0)

 check_arg_length(arg_value = eps,
                  arg_name = 'eps',
                  expected_length = 1)


 check_arg_type(arg_value = iter_max,
                arg_name = 'iter_max',
                expected_type = 'numeric')

 check_arg_is_integer(arg_value = iter_max,
                      arg_name = 'iter_max')

 check_arg_gteq(arg_value = iter_max,
                arg_name = 'iter_max',
                bound = 1)

 check_arg_length(arg_value = iter_max,
                  arg_name = 'iter_max',
                  expected_length = 1)


 check_arg_type(arg_value = pval_max,
                arg_name = 'pval_max',
                expected_type = 'numeric')

 check_arg_gt(arg_value = pval_max,
              arg_name = 'pval_max',
              bound = 0)

 check_arg_lteq(arg_value = pval_max,
                arg_name = 'pval_max',
                bound = 1)

 check_arg_length(arg_value = pval_max,
                  arg_name = 'pval_max',
                  expected_length = 1)


 check_arg_type(arg_value = do_scale,
                arg_name = 'do_scale',
                expected_type = 'logical')

 check_arg_length(arg_value = do_scale,
                  arg_name = 'do_scale',
                  expected_length = 1)

 if(!do_scale && iter_max > 1){
  stop("do_scale must be TRUE when iter_max > 1",
       call. = FALSE)
 }

}


#' Elastic net control
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

check_control_net <- function(alpha, df_target){

 check_arg_type(arg_value = alpha,
                arg_name = 'alpha',
                expected_type = 'numeric')

 check_arg_gteq(arg_value = alpha,
                arg_name = 'alpha',
                bound = 0)

 check_arg_lteq(arg_value = alpha,
                arg_name = 'alpha',
                bound = 1)

 check_arg_length(arg_value = alpha,
                  arg_name = 'alpha',
                  expected_length = 1)

 if(!is.null(df_target)){

  check_arg_type(arg_value = df_target,
                 arg_name = 'df_target',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = df_target,
                       arg_name = 'df_target')

 }

}
