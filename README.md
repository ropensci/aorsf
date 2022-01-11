
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
#>   0.80 0.2331220 0.1252519 0.04787173 0.3698876
#>    1.4 0.2496992 0.1448403 0.06142788 0.3871287
#>    3.5 0.3623194 0.2776568 0.15801605 0.5475111
#> 
#> -- age (VI Rank: 2) ----------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>     42 0.2704098 0.1417984 0.04608896 0.4440794
#>     50 0.2965021 0.1667501 0.05256073 0.5070508
#>     57 0.3273849 0.2118161 0.07401963 0.5635498
#> 
#> -- protime (VI Rank: 3) ------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>     10 0.2788376 0.1501081 0.05137325 0.4815079
#>     11 0.2910832 0.1592242 0.05324888 0.4945039
#>     11 0.3122619 0.1899187 0.06860109 0.5257119
#> 
#> -- copper (VI Rank: 4) -------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>     43 0.2615647 0.1451407 0.04761734 0.4448444
#>     74 0.2781534 0.1622355 0.05967552 0.4711515
#>    129 0.3295867 0.2258534 0.10952305 0.5317666
#> 
#> -- stage (VI Rank: 5) --------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>      1 0.2576518 0.1412697 0.04913651 0.4312446
#>      2 0.2671853 0.1405898 0.04736663 0.4673936
#>      3 0.2887433 0.1564222 0.05633898 0.5015073
#>      4 0.3355217 0.2158436 0.08858787 0.5588579
#> 
#>  Predicted risk at time t = 1788 for top 5 predictors
```
