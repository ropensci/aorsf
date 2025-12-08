# Oblique random forest control

Oblique random forest control

## Usage

``` r
orsf_control(
  tree_type,
  method,
  scale_x,
  ties,
  net_mix,
  target_df,
  max_iter,
  epsilon,
  ...
)

orsf_control_classification(
  method = "glm",
  scale_x = TRUE,
  net_mix = 0.5,
  target_df = NULL,
  max_iter = 20,
  epsilon = 1e-09,
  ...
)

orsf_control_regression(
  method = "glm",
  scale_x = TRUE,
  net_mix = 0.5,
  target_df = NULL,
  max_iter = 20,
  epsilon = 1e-09,
  ...
)

orsf_control_survival(
  method = "glm",
  scale_x = TRUE,
  ties = "efron",
  net_mix = 0.5,
  target_df = NULL,
  max_iter = 20,
  epsilon = 1e-09,
  ...
)
```

## Arguments

- tree_type:

  (*character*) the type of tree. Valid options are

  - "classification", i.e., categorical outcomes

  - "regression", i.e., continuous outcomes

  - "survival", i.e., time-to event outcomes

- method:

  (*character* or *function*) how to identify linear linear combinations
  of predictors. If `method` is a character value, it must be one of:

  - 'glm': linear, logistic, and cox regression

  - 'net': same as 'glm' but with penalty terms

  - 'random': random draw from uniform distribution

  If `method` is a *function*, it will be used to identify linear
  combinations of predictor variables. `method` must in this case accept
  three inputs named `x_node`, `y_node` and `w_node`, and should expect
  the following types and dimensions:

  - `x_node` (*matrix*; `n` rows, `p` columns)

  - `y_node` (*matrix*; `n` rows, `2` columns)

  - `w_node` (*matrix*; `n` rows, `1` column)

  In addition, `method` must return a matrix with p rows and 1 column.

- scale_x:

  (*logical*) if `TRUE`, values of predictors will be scaled prior to
  each instance of finding a linear combination of predictors, using
  summary values from the data in the current node of the decision tree.

- ties:

  (*character*) a character string specifying the method for tie
  handling. Only relevant when modeling survival outcomes and using a
  method that engages with tied outcome times. If there are no ties, all
  the methods are equivalent. Valid options are 'breslow' and 'efron'.
  The Efron approximation is the default because it is more accurate
  when dealing with tied event times and has similar computational
  efficiency compared to the Breslow method.

- net_mix:

  (*double*) The elastic net mixing parameter. A value of 1 gives the
  lasso penalty, and a value of 0 gives the ridge penalty. If multiple
  values of alpha are given, then a penalized model is fit using each
  alpha value prior to splitting a node.

- target_df:

  (*integer*) Preferred number of variables used in each linear
  combination. For example, with `mtry` of 5 and `target_df` of 3, we
  sample 5 predictors and look for the best linear combination using 3
  of them.

- max_iter:

  (*integer*) iteration continues until convergence (see `eps` above) or
  the number of attempted iterations is equal to `iter_max`.

- epsilon:

  (*double*) When using most modeling based method, iteration continues
  in the algorithm until the relative change in some kind of objective
  is less than `epsilon`, or the absolute change is less than
  `sqrt(epsilon)`.

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

an object of class `'orsf_control'`, which should be used as an input
for the `control` argument of
[orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md). Components
are:

- `tree_type`: type of trees to fit

- `lincomb_type`: method for linear combinations

- `lincomb_eps`: epsilon for convergence

- `lincomb_iter_max`: max iterations

- `lincomb_scale`: to scale or not.

- `lincomb_alpha`: mixing parameter

- `lincomb_df_target`: target degrees of freedom

- `lincomb_ties_method`: method for ties in survival time

- `lincomb_R_function`: R function for custom splits

## Details

Adjust `scale_x` *at your own risk*. Setting `scale_x = FALSE` will
reduce computation time but will also make the `orsf` model dependent on
the scale of your data, which is why the default value is `TRUE`.

## Examples

First we load some relevant packages

    set.seed(329730)
    suppressPackageStartupMessages({
     library(aorsf)
     library(survival)
     library(ranger)
     library(riskRegression)
    })

    ## Warning: package 'riskRegression' was built under R version 4.5.2

### Accelerated linear combinations

The accelerated ORSF ensemble is the default because it has a nice
balance of computational speed and prediction accuracy. It runs a single
iteration of Newton Raphson scoring on the Cox partial likelihood
function to find linear combinations of predictors.

    fit_accel <- orsf(pbc_orsf,
                      control = orsf_control_survival(),
                      formula = Surv(time, status) ~ . - id,
                      tree_seeds = 329)

### Linear combinations with Cox regression

Setting inputs in `orsf_control_survival` to scale the X matrix and
repeat iterations until convergence allows you to run Cox regression in
each non-terminal node of each survival tree, using the regression
coefficients to create linear combinations of predictors:

    control_cph <- orsf_control_survival(method = 'glm',
                                         scale_x = TRUE,
                                         max_iter = 20)

    fit_cph <- orsf(pbc_orsf,
                    control = control_cph,
                    formula = Surv(time, status) ~ . - id,
                    tree_seeds = 329)

### Linear combinations with penalized cox regression

Setting `method == 'net'` runs penalized Cox regression in each
non-terminal node of each survival tree. This can be really helpful if
you want to do feature selection within the node, but it is a lot slower
than the `'glm'` option.

    # select 3 predictors out of 5 to be used in
    # each linear combination of predictors.

    control_net <- orsf_control_survival(method = 'net', target_df = 3)

    fit_net <- orsf(pbc_orsf,
                    control = control_net,
                    formula = Surv(time, status) ~ . - id,
                    tree_seeds = 329)

### Linear combinations with your own function

In addition to the built-in methods, customized functions can be used to
identify linear combinations of predictors. We’ll demonstrate a few
here.

- The first uses random coefficients

    f_rando <- function(x_node, y_node, w_node){
     matrix(runif(ncol(x_node)), ncol=1)
    }

- The second derives coefficients from principal component analysis

    f_pca <- function(x_node, y_node, w_node) {

     # estimate two principal components.
     pca <- stats::prcomp(x_node, rank. = 2)
     # use the second principal component to split the node
     pca$rotation[, 1L, drop = FALSE]

    }

- The third uses `ranger()` inside of
  [`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md). This
  approach is very similar to a method known as reinforcement learning
  trees (see the `RLT` package), although our method of “muting” is very
  crude compared to the method proposed by Zhu et al. 

    f_rlt <- function(x_node, y_node, w_node){

     colnames(y_node) <- c('time', 'status')
     colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')

     data <- as.data.frame(cbind(y_node, x_node))

     if(nrow(data) <= 10)
      return(matrix(runif(ncol(x_node)), ncol = 1))

     fit <- ranger::ranger(data = data,
                           formula = Surv(time, status) ~ .,
                           num.trees = 25,
                           num.threads = 1,
                           min.node.size = 5,
                           importance = 'permutation')

     out <- sort(fit$variable.importance, decreasing = TRUE)

     # "mute" the least two important variables
     n_vars <- length(out)
     if(n_vars > 4){
       out[c(n_vars, n_vars-1)] <- 0
     }

     # ensure out has same variable order as input
     out <- out[colnames(x_node)]

     # protect yourself
     out[is.na(out)] <- 0

     matrix(out, ncol = 1)

    }

We can plug these functions into
[`orsf_control_custom()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_custom.md),
and then pass the result into
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md):

    fit_rando <- orsf(pbc_orsf,
                      Surv(time, status) ~ . - id,
                      control = orsf_control_survival(method = f_rando),
                      tree_seeds = 329)

    fit_pca <- orsf(pbc_orsf,
                    Surv(time, status) ~ . - id,
                    control = orsf_control_survival(method = f_pca),
                    tree_seeds = 329)

    fit_rlt <- orsf(pbc_orsf, time + status ~ . - id,
                    control = orsf_control_survival(method = f_rlt),
                    tree_seeds = 329)

So which fit seems to work best in this example? Let’s find out by
evaluating the out-of-bag survival predictions.

    risk_preds <- list(
     accel = fit_accel$pred_oobag,
     cph   = fit_cph$pred_oobag,
     net   = fit_net$pred_oobag,
     rando = fit_rando$pred_oobag,
     pca   = fit_pca$pred_oobag,
     rlt   = fit_rlt$pred_oobag
    )

    sc <- Score(object = risk_preds,
                formula = Surv(time, status) ~ 1,
                data = pbc_orsf,
                summary = 'IPA',
                times = fit_accel$pred_horizon)

The AUC values, from highest to lowest:

    sc$AUC$score[order(-AUC)]

    ##     model times       AUC         se     lower     upper
    ##    <fctr> <num>     <num>      <num>     <num>     <num>
    ## 1:    net  1788 0.9151649 0.02025057 0.8754745 0.9548553
    ## 2:  accel  1788 0.9095628 0.02143250 0.8675558 0.9515697
    ## 3:    cph  1788 0.9095628 0.02143250 0.8675558 0.9515697
    ## 4:    rlt  1788 0.9089871 0.02099354 0.8678406 0.9501337
    ## 5:  rando  1788 0.9062197 0.02148854 0.8641029 0.9483365
    ## 6:    pca  1788 0.8999479 0.02226683 0.8563057 0.9435901

And the indices of prediction accuracy:

    sc$Brier$score[order(-IPA), .(model, times, IPA)]

    ##         model times   IPA
    ##        <fctr> <num> <num>
    ## 1: Null model  1788    NA
    ## 2:      accel  1788    NA
    ## 3:        cph  1788    NA
    ## 4:        net  1788    NA
    ## 5:      rando  1788    NA
    ## 6:        pca  1788    NA
    ## 7:        rlt  1788    NA

From inspection,

- `net`, `accel`, and `rlt` have high discrimination and index of
  prediction accuracy.

- `rando` and `pca` do less well, but they aren’t bad.

## See also

linear combination control functions
[`orsf_control_cph()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_cph.md),
[`orsf_control_custom()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_custom.md),
[`orsf_control_fast()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_fast.md),
[`orsf_control_net()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_net.md)
