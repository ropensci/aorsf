
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
[![R-CMD-check](https://github.com/bcjaeger/aorsf/workflows/R-CMD-check/badge.svg)](https://github.com/bcjaeger/aorsf/actions)
<!-- badges: end -->

The goal of `aorsf` is to fit, interpret, and make predictions with
oblique random survival forests (ORSFs). The ‘a’ in the title stands for
accelerated. So why do ORSFs need to be accelerated? Oblique decision
trees are often more accurate but slower to fit compared to their
axis-based counterparts. The issue of higher computation time is
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
#> -- age (VI Rank: 1) ----------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>   41.5 0.2706164 0.1337351 0.04630290 0.4612794
#>   49.7 0.2975084 0.1661233 0.05237698 0.5170854
#>   56.6 0.3284048 0.2116927 0.07674439 0.5457127
#> 
#> -- bili (VI Rank: 2) ---------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>   0.80 0.2320563 0.1234708 0.04876897 0.3530933
#>   1.40 0.2504298 0.1421328 0.06205604 0.3746472
#>   3.52 0.3650838 0.2800719 0.16607641 0.5189161
#> 
#> -- sex (VI Rank: 3) ----------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>      m 0.3536771 0.2485654 0.11330040 0.5797451
#>      f 0.2928458 0.1494558 0.05061492 0.5113725
#> 
#> -- spiders (VI Rank: 4) ------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>      0 0.2881293 0.1477039 0.04940008 0.4931721
#>      1 0.3327476 0.2142004 0.08124437 0.5465791
#> 
#> -- copper (VI Rank: 5) -------------------------
#> 
#>        |---------------- risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>   42.8 0.2619422 0.1423453 0.04838753 0.4504369
#>   74.0 0.2789013 0.1513143 0.05669541 0.4814084
#>    129 0.3313322 0.2132656 0.10764336 0.5389831
#> 
#>  Predicted risk at time t = 1788 for top 5 predictors
```

## References

Byron C. Jaeger, D. Leann Long, Dustin M. Long, Mario Sims, Jeff M.
Szychowski, Yuan-I Min, Leslie A. Mcclure, George Howard, Noah Simon
(2019). Oblique Random Survival Forests. Ann. Appl. Stat. 13(3):
1847-1883. URL <https://doi.org/10.1214/19-AOAS1261> DOI:
10.1214/19-AOAS1261
