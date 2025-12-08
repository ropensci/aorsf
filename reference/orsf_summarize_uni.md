# Univariate summary

Summarize the univariate information from an ORSF object

## Usage

``` r
orsf_summarize_uni(
  object,
  n_variables = NULL,
  pred_horizon = NULL,
  pred_type = NULL,
  importance = NULL,
  class = NULL,
  verbose_progress = FALSE,
  ...
)
```

## Arguments

- object:

  (*ObliqueForest*) a trained oblique random forest object (see
  [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md)).

- n_variables:

  (*integer*) how many variables should be summarized? Setting this
  input to a lower number will reduce computation time.

- pred_horizon:

  (*double*) Only relevent for survival forests. A value or vector
  indicating the time(s) that predictions will be calibrated to. E.g.,
  if you were predicting risk of incident heart failure within the next
  10 years, then `pred_horizon = 10`. `pred_horizon` can be `NULL` if
  `pred_type` is `'mort'`, since mortality predictions are aggregated
  over all event times

- pred_type:

  (*character*) the type of predictions to compute. Valid Valid options
  for survival are:

  - 'risk' : probability of having an event at or before `pred_horizon`.

  - 'surv' : 1 - risk.

  - 'chf': cumulative hazard function

  - 'mort': mortality prediction

  - 'time': survival time prediction

  For classification:

  - 'prob': probability for each class

  For regression:

  - 'mean': predicted mean, i.e., the expected value

- importance:

  (*character*) Indicate method for variable importance:

  - 'none': no variable importance is computed.

  - 'anova': compute analysis of variance (ANOVA) importance

  - 'negate': compute negation importance

  - 'permute': compute permutation importance

- class:

  (*character*) only relevant for classification forests. If `NULL` (the
  default), summary statistics are returned for all classes in the
  outcome, and printed summaries will show the last class in the class
  levels. To specify a single class to summarize, indicate the name of
  the class with `class`. E.g., if the categorical outcome has class
  levels A, B, and C, then using `class = "A"` will restrict output to
  class A.

  For details on these methods, see
  [orsf_vi](https://bcjaeger.github.io/aorsf/reference/orsf_vi.md).

- verbose_progress:

  (*logical*) if `TRUE`, progress will be printed to console. If `FALSE`
  (the default), nothing will be printed.

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

an object of class 'orsf_summary', which includes data on

- importance of individual predictors.

- expected values of predictions at specific values of predictors.

## Details

If `pred_horizon` is left unspecified, the median value of the
time-to-event variable in `object`'s training data will be used. It is
recommended to always specify your own prediction horizon, as the median
time may not be an especially meaningful horizon to compute predicted
risk values at.

If `object` already has variable importance values, you can safely
bypass the computation of variable importance in this function by
setting importance = 'none'.

## See also

as.data.table.orsf_summary_uni

## Examples

``` r
object <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 25)

# since anova importance was used to make object, it is also
# used for ranking variables in the summary, unless we specify
# a different type of importance

orsf_summarize_uni(object, n_variables = 2)
#> 
#> -- ascites (VI Rank: 1) -------------------------
#> 
#>         |---------------- Risk ----------------|
#>   Value      Mean    Median     25th %    75th %
#>  <char>     <num>     <num>      <num>     <num>
#>       0 0.3111992 0.1938980 0.05124179 0.5400958
#>       1 0.4670078 0.4462577 0.26074220 0.6430327
#> 
#> -- bili (VI Rank: 2) ----------------------------
#> 
#>         |---------------- Risk ----------------|
#>   Value      Mean    Median     25th %    75th %
#>  <char>     <num>     <num>      <num>     <num>
#>    0.60 0.2302508 0.1432700 0.03912098 0.3734274
#>    0.80 0.2372870 0.1481850 0.04209576 0.3784865
#>    1.40 0.2618433 0.1830631 0.06358040 0.4004983
#>    3.52 0.4047835 0.3604112 0.22291673 0.5523087
#>    7.25 0.5121593 0.4893694 0.34706191 0.6712866
#> 
#>  Predicted risk at time t = 1788 for top 2 predictors 

# if we want to summarize object according to variables
# ranked by negation importance, we can compute negation
# importance within orsf_summarize_uni() as follows:

orsf_summarize_uni(object, n_variables = 2, importance = 'negate')
#> 
#> -- bili (VI Rank: 1) ----------------------------
#> 
#>         |---------------- Risk ----------------|
#>   Value      Mean    Median     25th %    75th %
#>  <char>     <num>     <num>      <num>     <num>
#>    0.60 0.2302508 0.1432700 0.03912098 0.3734274
#>    0.80 0.2372870 0.1481850 0.04209576 0.3784865
#>    1.40 0.2618433 0.1830631 0.06358040 0.4004983
#>    3.52 0.4047835 0.3604112 0.22291673 0.5523087
#>    7.25 0.5121593 0.4893694 0.34706191 0.6712866
#> 
#> -- sex (VI Rank: 2) -----------------------------
#> 
#>         |---------------- Risk ----------------|
#>   Value      Mean    Median     25th %    75th %
#>  <char>     <num>     <num>      <num>     <num>
#>       m 0.3429356 0.2383294 0.10642774 0.5506692
#>       f 0.3142440 0.1841889 0.05244527 0.5475021
#> 
#>  Predicted risk at time t = 1788 for top 2 predictors 

# for multi-category fits, you can specify which class
# you want to summarize:

object =  orsf(species ~ ., data = penguins_orsf, n_tree = 25)

orsf_summarize_uni(object, class = "Adelie", n_variables = 1)
#> 
#> -- bill_length_mm (VI Rank: 1) ------------------
#> 
#>         |------------- Probability -------------|
#>   Value      Mean    Median     25th %    75th %
#>  <char>     <num>     <num>      <num>     <num>
#>    36.6 0.6919210 0.8484422 0.37763690 0.9664300
#>    39.5 0.6614623 0.7927307 0.33801079 0.9613584
#>    44.5 0.4204350 0.4635360 0.06783333 0.6781295
#>    48.6 0.2011209 0.1551832 0.01302446 0.3397603
#>    50.8 0.1701162 0.1212186 0.01023102 0.2895517
#> 
#>  Predicted probability for top 1 predictors 

```
