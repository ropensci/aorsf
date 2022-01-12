
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
decision trees are often more accurate but slower to fit compared to
their axis-based counterparts. The issue of higher computation time is
compounded for survival decision trees, which usually require more
computing than classification or regression trees. To make ORSF more
accessible and able to engage with larger datasets, `aorsf` applies
strategies to cut down computing time without sacrificing prediction
accuracy.

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
#>   0.80 0.2303045 0.1159666 0.04540949 0.3738375
#>   1.40 0.2501020 0.1314716 0.06338070 0.3913647
#>   3.52 0.3614245 0.2718542 0.15921464 0.5318559
#> 
#> -- age (VI Rank: 2) ----------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>   41.5 0.2700114 0.1346919 0.04416058 0.4529193
#>   49.7 0.2957704 0.1642519 0.05406492 0.4954041
#>   56.6 0.3254694 0.2082042 0.07406272 0.5352754
#> 
#> -- ascites (VI Rank: 3) ------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>      0 0.2917184 0.1545333 0.05076274 0.5028834
#>      1 0.4627266 0.3794409 0.25965444 0.6429872
#> 
#> -- protime (VI Rank: 4) ------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>   10.0 0.2787354 0.1502139 0.04815788 0.4908722
#>   10.6 0.2916406 0.1630598 0.05422335 0.5144108
#>   11.2 0.3119439 0.1838502 0.07078713 0.5324805
#> 
#> -- albumin (VI Rank: 5) ------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>   3.31 0.3128104 0.1812829 0.05640159 0.5564690
#>   3.54 0.2890683 0.1506630 0.04749213 0.5105134
#>   3.77 0.2750514 0.1459699 0.04861506 0.4788465
#> 
#>  Predicted risk at time t = 1788 for top 5 predictors
```

## References

Byron C. Jaeger, D. Leann Long, Dustin M. Long, Mario Sims, Jeff M.
Szychowski, Yuan-I Min, Leslie A. Mcclure, George Howard, Noah Simon
(2019). Oblique Random Survival Forests. Ann. Appl. Stat. 13(3):
1847-1883. URL <https://doi.org/10.1214/19-AOAS1261> DOI:
10.1214/19-AOAS1261
