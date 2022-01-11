
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
[![R-CMD-check](https://github.com/bcjaeger/aorsf/workflows/R-CMD-check/badge.svg)](https://github.com/bcjaeger/aorsf/actions)
<!-- badges: end -->

The goal of `aorsf` is to fit, interpret, and make predictions with
oblique random survival forests (ORSFs).

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
values of those predictors.

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
#>   0.80 0.2325340 0.1265260 0.04848649 0.3757404
#>    1.4 0.2506499 0.1408313 0.06174594 0.3890790
#>    3.5 0.3672571 0.2771605 0.16780738 0.5592806
#> 
#> -- age (VI Rank: 2) ----------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>     42 0.2742214 0.1506950 0.04545622 0.4585691
#>     50 0.2977520 0.1759675 0.05254730 0.5043843
#>     57 0.3268951 0.2143951 0.07200863 0.5633268
#> 
#> -- ascites (VI Rank: 3) ------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>      0 0.2953194 0.1688652 0.05285233 0.5296210
#>      1 0.4711030 0.3923678 0.28036067 0.6666691
#> 
#> -- protime (VI Rank: 4) ------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>     10 0.2802406 0.1542821 0.04971174 0.5090750
#>     11 0.2932556 0.1708379 0.05749268 0.5189746
#>     11 0.3170823 0.1904355 0.07086467 0.5495897
#> 
#> -- spiders (VI Rank: 5) ------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>      0 0.2903928 0.1590619 0.04953802 0.5106214
#>      1 0.3361667 0.2161500 0.09189305 0.5474268
#> 
#>  Predicted risk at time t = 1788 for top 5 predictors
```

The term ‘uni’ is short for univariate.
