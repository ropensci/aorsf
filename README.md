
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
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
#>   0.80 0.2334629 0.1250566 0.04667434 0.3578025
#>    1.4 0.2511046 0.1389551 0.05616435 0.3905442
#>    3.5 0.3692912 0.2931586 0.16049991 0.5355958
#> 
#> -- copper (VI Rank: 2) -------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>     43 0.2636023 0.1513536 0.04835872 0.4403300
#>     74 0.2785638 0.1682814 0.05766171 0.4733305
#>    129 0.3316113 0.2225723 0.10986226 0.5310440
#> 
#> -- age (VI Rank: 3) ----------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>     42 0.2704019 0.1376899 0.04833119 0.4560358
#>     50 0.2962462 0.1653618 0.05543100 0.5055893
#>     57 0.3272574 0.2089887 0.07685078 0.5481599
#> 
#> -- ast (VI Rank: 4) ----------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>     82 0.2835042 0.1565048 0.04901267 0.5133886
#>    117 0.2976624 0.1643752 0.05407187 0.5352556
#>    153 0.3194210 0.1824407 0.07580561 0.5691429
#> 
#> -- ascites (VI Rank: 5) ------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>      0 0.2941116 0.1637565 0.05421877 0.5249933
#>      1 0.4732697 0.3991843 0.28792132 0.6352829
#> 
#>  Predicted risk at time t = 1788 for top 5 predictors
```

The term ‘uni’ is short for univariate.
