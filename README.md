
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
[![R-CMD-check](https://github.com/bcjaeger/aorsf/workflows/R-CMD-check/badge.svg)](https://github.com/bcjaeger/aorsf/actions)
[![pkgcheck](https://github.com/bcjaeger/aorsf/workflows/pkgcheck/badge.svg)](https://github.com/bcjaeger/aorsf/actions?query=workflow%3Apkgcheck)
[![Status at rOpenSci Software Peer
Review](https://badges.ropensci.org/532_status.svg)](https://github.com/ropensci/software-review/issues/532)
<!-- badges: end -->

`aorsf` provides optimized software to fit, interpret, and make
predictions with oblique random survival forests (ORSFs).

## Why aorsf?

-   over 400 times faster than `obliqueRSF`.

-   accurate predictions for time-to-event outcomes.

-   negation importance, a novel technique to estimate variable
    importance for ORSFs.

-   intuitive API with formula based interface.

-   extensive input checks + informative error messages.

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
#>           OOB stat value: 0.84
#>            OOB stat type: Harrell's C-statistic
#> 
#> -----------------------------------------
```

How about interpreting the fit?

-   use `orsf_vi_negate()` and `orsf_vi_anova()` for variable importance

    ``` r
    orsf_vi_negate(fit)
    #>          bili           age        copper       ascites       protime 
    #>  1.390915e-02  1.213795e-02  8.387164e-03  6.511773e-03  6.199208e-03 
    #>         stage           sex        hepato         edema       albumin 
    #>  4.428006e-03  4.375912e-03  3.490310e-03  3.170303e-03  2.240050e-03 
    #>       spiders           trt      alk.phos           ast          trig 
    #>  1.041884e-03  1.562826e-04  5.209419e-05 -1.198166e-03 -1.250260e-03 
    #>      platelet          chol 
    #> -1.406543e-03 -5.782455e-03
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>     bili      mean        lwr      medn       upr
    #>    <int>     <num>      <num>     <num>     <num>
    #> 1:     1 0.2327332 0.01715402 0.1278394 0.8629424
    #> 2:     2 0.2803282 0.03730412 0.1726425 0.8873881
    #> 3:     3 0.3399960 0.06581874 0.2534968 0.9011461
    #> 4:     4 0.3936117 0.10559414 0.3220519 0.9190984
    #> 5:     5 0.4360391 0.13796185 0.3741413 0.9339884
    ```

-   use `orsf_summarize_uni()` to show the top predictor variables in an
    ORSF model and the expected predicted risk at specific values of
    those predictors. (The term ‘uni’ is short for univariate.)

    ``` r
    # take a look at the top 5 variables 
    # for continuous predictors, see expected risk at 25/50/75th quantile
    # for categorical predictors, see expected risk in specified category

    orsf_summarize_uni(object = fit, n_variables = 5)
    #> 
    #> -- bili (VI Rank: 1) ----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>    0.80 0.2276726 0.1248961 0.05111495 0.3557978
    #>    1.40 0.2473343 0.1394002 0.06143315 0.3806993
    #>    3.52 0.3692963 0.2953896 0.17716271 0.5528412
    #> 
    #> -- age (VI Rank: 2) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>    41.5 0.2713199 0.1407412 0.04671284 0.4618216
    #>    49.7 0.2984502 0.1723371 0.05323084 0.5051649
    #>    56.6 0.3316576 0.2264977 0.07716886 0.5633089
    #> 
    #> -- copper (VI Rank: 3) --------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>    42.8 0.2640370 0.1344469 0.04820217 0.4572128
    #>    74.0 0.2812038 0.1607180 0.05813276 0.4765200
    #>     129 0.3299795 0.2238479 0.09849015 0.5367875
    #> 
    #> -- ascites (VI Rank: 4) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean   Median    25th %    75th %
    #>  <char>     <num>    <num>     <num>     <num>
    #>       0 0.2942636 0.160008 0.0501474 0.5045309
    #>       1 0.4547808 0.373920 0.2533511 0.6365145
    #> 
    #> -- protime (VI Rank: 5) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>    10.0 0.2820317 0.1483283 0.04867649 0.4915174
    #>    10.6 0.2947320 0.1590040 0.05496065 0.5056502
    #>    11.2 0.3154625 0.1938762 0.06724748 0.5468626
    #> 
    #>  Predicted risk at time t = 1788 for top 5 predictors
    ```

## References

Byron C. Jaeger, D. Leann Long, Dustin M. Long, Mario Sims, Jeff M.
Szychowski, Yuan-I Min, Leslie A. Mcclure, George Howard, Noah Simon
(2019). Oblique Random Survival Forests. Ann. Appl. Stat. 13(3):
1847-1883. URL <https://doi.org/10.1214/19-AOAS1261> DOI:
10.1214/19-AOAS1261

## Funding

This software receives financial support from the Center for Biomedical
Informatics, Wake Forest School of Medicine.
