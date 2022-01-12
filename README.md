
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
[![R-CMD-check](https://github.com/bcjaeger/aorsf/workflows/R-CMD-check/badge.svg)](https://github.com/bcjaeger/aorsf/actions)
<!-- badges: end -->

The goal of `aorsf` is to fit, interpret, and make predictions with
oblique random survival forests (ORSFs). The ‘a’ in the title of `aorsf`
stands for accelerated. So why do ORSFs need to be accelerated? Oblique
decision trees are notoriously slow, and the issue of higher computation
time is compounded by the fast that survival trees may also be slow. To
make ORSF more accessible and able to engage with larger datasets,
`aorsf` applies strategies to cut down computing time without
sacrificing prediction accuracy.

## Installation

You can install the development version of aorsf from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("bcjaeger/aorsf")
```

## Example

The `orsf()` function is used to fit ORSFs. Printing the output from
`orsf()` will give some descriptive statistics about the ensemble.

``` r
library(aorsf)

fit <- orsf(data_train = pbc_orsf,
            formula = Surv(time, status) ~ . - id)

print(fit)
#> ---------- Oblique random survival forest
#> 
#>           N observations: 276
#>                 N events: 111
#>                  N trees: 500
#>       N predictors total: 17
#>    N predictors per node: 5
#>  Average leaves per tree: 19
#> Min observations in leaf: 5
#>       Min events in leaf: 1
#>          OOB C-statistic: 0.84
#> 
#> -----------------------------------------
```

How about interpreting the fit? There are several functions for this:
`orsf_vi()` for variable importance, `orsf_interaction()` for two-way
variable interactions, and `orsf_pd_ice()` or `orsf_pd_summary()` for
individual or aggregated partial dependence values. However,
`orsf_summarize_uni()` is the most convenient way to assess top
predictor variables in ORSF and the expected predicted risk at specific
values of those predictors. The term ‘uni’ is short for univariate.

``` r
# take a look at the top 5 variables 
# for continuous predictors, see expected risk at 25/50/75th quantile
# for categorical predictors, see expected risk in each category

orsf_summarize_uni(object = fit, n_variables = 5)
#> 
#> -- bili (VI Rank: 1) ---------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>   0.80 0.2327223 0.1187105 0.04920697 0.3532136
#>    1.4 0.2503352 0.1355363 0.06129695 0.3794809
#>    3.5 0.3617920 0.2801387 0.16335693 0.5145361
#> 
#> -- age (VI Rank: 2) ----------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>     42 0.2719462 0.1368686 0.04647107 0.4459963
#>     50 0.2955702 0.1644157 0.05180857 0.5115135
#>     57 0.3237065 0.2026756 0.07015707 0.5451979
#> 
#> -- ascites (VI Rank: 3) ------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>      0 0.2914747 0.1571831 0.05285387 0.5241465
#>      1 0.4699443 0.3873810 0.28113338 0.6587456
#> 
#> -- stage (VI Rank: 4) --------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>      1 0.2579344 0.1377085 0.04771583 0.4422156
#>      2 0.2657589 0.1446005 0.04661316 0.4628640
#>      3 0.2884751 0.1683806 0.05461591 0.5014843
#>      4 0.3347969 0.2075534 0.08498729 0.5744516
#> 
#> -- copper (VI Rank: 5) -------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>     43 0.2599877 0.1425730 0.05012715 0.4470400
#>     74 0.2765080 0.1562001 0.05850182 0.4536671
#>    129 0.3319711 0.2266852 0.11518515 0.5382018
#> 
#>  Predicted risk at time t = 1788 for top 5 predictors
```
