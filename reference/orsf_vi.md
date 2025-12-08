# Variable Importance

Estimate the importance of individual predictor variables using oblique
random forests.

## Usage

``` r
orsf_vi(
  object,
  group_factors = TRUE,
  importance = NULL,
  oobag_fun = NULL,
  n_thread = NULL,
  verbose_progress = NULL,
  ...
)

orsf_vi_negate(
  object,
  group_factors = TRUE,
  oobag_fun = NULL,
  n_thread = NULL,
  verbose_progress = NULL,
  ...
)

orsf_vi_permute(
  object,
  group_factors = TRUE,
  oobag_fun = NULL,
  n_thread = NULL,
  verbose_progress = NULL,
  ...
)

orsf_vi_anova(object, group_factors = TRUE, verbose_progress = NULL, ...)
```

## Arguments

- object:

  (*ObliqueForest*) a trained oblique random forest object (see
  [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md)).

- group_factors:

  (*logical*) if `TRUE`, the importance of factor variables will be
  reported overall by aggregating the importance of individual levels of
  the factor. If `FALSE`, the importance of individual factor levels
  will be returned.

- importance:

  (*character*) Indicate method for variable importance:

  - 'anova': compute analysis of variance (ANOVA) importance

  - 'negate': compute negation importance

  - 'permute': compute permutation importance

- oobag_fun:

  (*function*) to be used for evaluating out-of-bag prediction accuracy
  after negating coefficients (if importance = 'negate') or permuting
  the values of a predictor (if importance = 'permute')

  - When `oobag_fun = NULL` (the default), the evaluation statistic is
    selected based on tree type

  - survival: Harrell's C-statistic (1982)

  - classification: Area underneath the ROC curve (AUC-ROC)

  - regression: Traditional prediction R-squared

  - if you use your own `oobag_fun` note the following:

    - `oobag_fun` should have three inputs: `y_mat`, `w_vec`, and
      `s_vec`

    - For survival trees, `y_mat` should be a two column matrix with
      first column named 'time' and second named 'status'. For
      classification trees, `y_mat` should be a matrix with number of
      columns = number of distinct classes in the outcome. For
      regression, `y_mat` should be a matrix with one column.

    - `s_vec` is a numeric vector containing predictions

    - `oobag_fun` should return a numeric output of length 1

    - the same `oobag_fun` should have been used when you created
      `object` so that the initial value of out-of-bag prediction
      accuracy is consistent with the values that will be computed while
      variable importance is estimated.

  For more details, see the out-of-bag
  [vignette](https://docs.ropensci.org/aorsf/articles/oobag.html).

- n_thread:

  (*integer*) number of threads to use while computing predictions.
  Default is 0, which allows a suitable number of threads to be used
  based on availability.

- verbose_progress:

  (*logical*) if `TRUE`, progress messages are printed in the console.
  If `FALSE` (the default), nothing is printed.

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

`orsf_vi` functions return a named numeric vector.

- Names of the vector are the predictor variables used by `object`

- Values of the vector are the estimated importance of the given
  predictor.

The returned vector is sorted from highest to lowest value, with higher
values indicating higher importance.

## Details

When an `ObliqueForest` object is grown with importance = 'anova',
'negate', or 'permute', the output will have a vector of importance
values based on the requested type of importance. However, `orsf_vi()`
can be used to compute variable importance after growing a forest or to
compute a different type of importance.

`orsf_vi()` is a general purpose function to extract or compute variable
importance estimates from an `ObliqueForest` object (see
[orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md)).
`orsf_vi_negate()`, `orsf_vi_permute()`, and `orsf_vi_anova()` are
wrappers for `orsf_vi()`. The way these functions work depends on
whether the `object` they are given already has variable importance
estimates in it or not (see examples).

## Variable importance methods

**negation importance**: Each variable is assessed separately by
multiplying the variable's coefficients by -1 and then determining how
much the model's performance changes. The worse the model's performance
after negating coefficients for a given variable, the more important the
variable. This technique is promising b/c it does not require
permutation and it emphasizes variables with larger coefficients in
linear combinations, but it is also relatively new and hasn't been
studied as much as permutation importance. See Jaeger, (2023) for more
details on this technique.

**permutation importance**: Each variable is assessed separately by
randomly permuting the variable's values and then determining how much
the model's performance changes. The worse the model's performance after
permuting the values of a given variable, the more important the
variable. This technique is flexible, intuitive, and frequently used. It
also has several [known
limitations](https://christophm.github.io/interpretable-ml-book/feature-importance.html#disadvantages-9)

**analysis of variance (ANOVA) importance**: A p-value is computed for
each coefficient in each linear combination of variables in each
decision tree. Importance for an individual predictor variable is the
proportion of times a p-value for its coefficient is \< 0.01. This
technique is very efficient computationally, but may not be as effective
as permutation or negation in terms of selecting signal over noise
variables. See [Menze,
2011](https://link.springer.com/chapter/10.1007/978-3-642-23783-6_29)
for more details on this technique.

## Examples

### ANOVA importance

The default variable importance technique, ANOVA, is calculated while
you fit an oblique random forest ensemble.

    fit <- orsf(pbc_orsf, Surv(time, status) ~ . - id)

    fit

    ## ---------- Oblique random survival forest
    ##
    ##      Linear combinations: Accelerated Cox regression
    ##           N observations: 276
    ##                 N events: 111
    ##                  N trees: 500
    ##       N predictors total: 17
    ##    N predictors per node: 5
    ##  Average leaves per tree: 21.022
    ## Min observations in leaf: 5
    ##       Min events in leaf: 1
    ##           OOB stat value: 0.84
    ##            OOB stat type: Harrell's C-index
    ##      Variable importance: anova
    ##
    ## -----------------------------------------

ANOVA is the default because it is fast, but it may not be as decisive
as the permutation and negation techniques for variable selection.

### Raw VI values

the ‘raw’ variable importance values can be accessed from the fit object

    fit$get_importance_raw()

    ##                   [,1]
    ## trt_placebo 0.06355042
    ## age         0.23259259
    ## sex_f       0.14700432
    ## ascites_1   0.46791708
    ## hepato_1    0.14349776
    ## spiders_1   0.17371938
    ## edema_0.5   0.17459191
    ## edema_1     0.51197605
    ## bili        0.40590758
    ## chol        0.17666667
    ## albumin     0.25972156
    ## copper      0.28840580
    ## alk.phos    0.10614251
    ## ast         0.18327491
    ## trig        0.12815626
    ## platelet    0.09265648
    ## protime     0.22656250
    ## stage       0.20264766

these are ‘raw’ because values for factors have not been aggregated into
a single value. Currently there is one value for k-1 levels of a k level
factor. For example, you can see edema_1 and edema_0.5 in the importance
values above because edema is a factor variable with levels of 0, 0.5,
and 1.

### Collapse VI across factor levels

To get aggregated values across all levels of each factor,

- access the `importance` element from the `orsf` fit:

      # this assumes you used group_factors = TRUE in orsf()
      fit$importance

      ##    ascites       bili      edema     copper    albumin        age    protime
      ## 0.46791708 0.40590758 0.31115216 0.28840580 0.25972156 0.23259259 0.22656250
      ##      stage        ast       chol    spiders        sex     hepato       trig
      ## 0.20264766 0.18327491 0.17666667 0.17371938 0.14700432 0.14349776 0.12815626
      ##   alk.phos   platelet        trt
      ## 0.10614251 0.09265648 0.06355042

- use `orsf_vi()` with group_factors set to `TRUE` (the default)

      orsf_vi(fit)

      ##    ascites       bili      edema     copper    albumin        age    protime
      ## 0.46791708 0.40590758 0.31115216 0.28840580 0.25972156 0.23259259 0.22656250
      ##      stage        ast       chol    spiders        sex     hepato       trig
      ## 0.20264766 0.18327491 0.17666667 0.17371938 0.14700432 0.14349776 0.12815626
      ##   alk.phos   platelet        trt
      ## 0.10614251 0.09265648 0.06355042

Note that you can make the default returned importance values ungrouped
by setting `group_factors` to `FALSE` in the `orsf_vi` functions or the
`orsf` function.

### Add VI to an oblique random forest

You can fit an oblique random forest without VI, then add VI later

    fit_no_vi <- orsf(pbc_orsf,
                      Surv(time, status) ~ . - id,
                      importance = 'none')

    # Note: you can't call orsf_vi_anova() on fit_no_vi because anova
    # VI can only be computed while the forest is being grown.

    orsf_vi_negate(fit_no_vi)

    ##        bili      copper         sex     protime         age       stage
    ## 0.130439814 0.051880867 0.038308025 0.025115249 0.023826061 0.020354822
    ##     albumin     ascites        chol         ast     spiders      hepato
    ## 0.019997729 0.015918292 0.013320469 0.010086726 0.007409116 0.007326714
    ##       edema         trt    alk.phos        trig    platelet
    ## 0.006844435 0.003214544 0.002517057 0.002469545 0.001056829

    orsf_vi_permute(fit_no_vi)

    ##          bili        copper           age       ascites       protime
    ##  0.0592069141  0.0237362075  0.0136479213  0.0130805894  0.0123091354
    ##         stage       albumin          chol        hepato           ast
    ##  0.0117177661  0.0106414724  0.0064501213  0.0058813969  0.0057753740
    ##         edema       spiders           sex          trig      platelet
    ##  0.0052171180  0.0048427005  0.0023386947  0.0017883700  0.0013533691
    ##      alk.phos           trt
    ##  0.0006492029 -0.0009921507

### Oblique random forest and VI all at once

fit an oblique random forest and compute vi at the same time

    fit_permute_vi <- orsf(pbc_orsf,
                           Surv(time, status) ~ . - id,
                           importance = 'permute')

    # get the vi instantly (i.e., it doesn't need to be computed again)
    orsf_vi_permute(fit_permute_vi)

    ##          bili        copper       ascites       protime       albumin
    ##  0.0571305446  0.0243657146  0.0138318057  0.0133401675  0.0130746154
    ##           age         stage          chol           ast       spiders
    ##  0.0123610374  0.0102963203  0.0077895394  0.0075250059  0.0048628813
    ##         edema        hepato           sex      platelet          trig
    ##  0.0046003168  0.0039818730  0.0016891584  0.0012767063  0.0007324402
    ##      alk.phos           trt
    ##  0.0005128897 -0.0014443967

You can still get negation VI from this fit, but it needs to be computed

    orsf_vi_negate(fit_permute_vi)

    ##        bili      copper         sex     protime       stage         age
    ## 0.123331760 0.052544318 0.037291358 0.024977898 0.023239189 0.021934511
    ##     albumin     ascites        chol         ast     spiders       edema
    ## 0.020586632 0.014229536 0.014053040 0.012227048 0.007643156 0.006832766
    ##      hepato         trt    alk.phos        trig    platelet
    ## 0.006301693 0.004348705 0.002371797 0.002309396 0.001347035

### Custom functions for VI

The default prediction accuracy functions work well most of the time:

    fit_standard <- orsf(penguins_orsf, bill_length_mm ~ ., tree_seeds = 1)

    # Default method for prediction accuracy with VI is R-squared
    orsf_vi_permute(fit_standard)

    ##           species flipper_length_mm       body_mass_g     bill_depth_mm
    ##      0.3725898166      0.3261834607      0.2225730676      0.1026569498
    ##            island               sex              year
    ##      0.0876071687      0.0844807334      0.0006978493

But sometimes you want to do something specific and the defaults just
won’t work. For these cases, you can compute VI with any function you’d
like to measure prediction accuracy by supplying a valid function to the
`oobag_fun` input. For example, we use mean absolute error below. Higher
values are considered good when `aorsf` computes prediction accuracy, so
we make our function return a pseudo R-squared based on mean absolute
error:

    rsq_mae <- function(y_mat, w_vec, s_vec){

     mae_standard <- mean(abs((y_mat - mean(y_mat)) * w_vec))
     mae_fit <- mean(abs((y_mat - s_vec) * w_vec))

     1 - mae_fit / mae_standard

    }

    fit_custom <- orsf_update(fit_standard, oobag_fun = rsq_mae)

    # not much changes, but the difference between variables shrinks
    # and the ordering of sex and island has swapped
    orsf_vi_permute(fit_custom)

    ##           species flipper_length_mm       body_mass_g     bill_depth_mm
    ##       0.206951751       0.193248912       0.140899603       0.076759148
    ##               sex            island              year
    ##       0.073042331       0.050851073       0.003633365

## References

1.  Harrell, E F, Califf, M R, Pryor, B D, Lee, L K, Rosati, A R (1982).
    "Evaluating the yield of medical tests." *Jama*, *247*(18),
    2543-2546.

2.  Breiman, Leo (2001). "Random Forests." *Machine Learning*, *45*(1),
    5-32. ISSN 1573-0565.

3.  Menze, H B, Kelm, Michael B, Splitthoff, N D, Koethe, Ullrich,
    Hamprecht, A F (2011). "On oblique random forests." In *Machine
    Learning and Knowledge Discovery in Databases: European Conference,
    ECML PKDD 2011, Athens, Greece, September 5-9, 2011, Proceedings,
    Part II 22*, 453-469. Springer.

4.  Jaeger BC, Welden S, Lenoir K, Speiser JL, Segar MW, Pandey A,
    Pajewski NM (2023). "Accelerated and interpretable oblique random
    survival forests." *Journal of Computational and Graphical
    Statistics*, 1-16.
