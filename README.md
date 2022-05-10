
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

-   over 500 times faster than `obliqueRSF`.

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
#> Warning: package 'aorsf' was built under R version 4.1.3

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
    #>          bili           age           sex        copper       protime 
    #>  0.0169827047  0.0141175245  0.0087518233  0.0075536570  0.0065638675 
    #>       ascites         stage       spiders       albumin        hepato 
    #>  0.0055740779  0.0049489477  0.0042196291  0.0036465930  0.0032819337 
    #>         edema           ast          trig          chol      alk.phos 
    #>  0.0029582056  0.0022921442  0.0002083767 -0.0002604709 -0.0020837675 
    #>           trt      platelet 
    #> -0.0021879558 -0.0024484268
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>     bili      mean        lwr      medn       upr
    #>    <int>     <num>      <num>     <num>     <num>
    #> 1:     1 0.2372611 0.01574818 0.1294677 0.8765468
    #> 2:     2 0.2823997 0.03808830 0.1751359 0.9008787
    #> 3:     3 0.3342751 0.06253403 0.2385166 0.9125669
    #> 4:     4 0.3846965 0.08521571 0.3106505 0.9235278
    #> 5:     5 0.4281614 0.11914911 0.3608596 0.9266794
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
    #>    0.80 0.2332029 0.1229387 0.04171414 0.3596673
    #>    1.40 0.2503437 0.1396634 0.05286906 0.3873875
    #>    3.52 0.3602802 0.2742367 0.14994833 0.5389175
    #> 
    #> -- age (VI Rank: 2) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>    41.5 0.2721148 0.1360468 0.04008261 0.4607133
    #>    49.7 0.2969165 0.1610665 0.05026247 0.5187667
    #>    56.6 0.3266064 0.2078738 0.07281189 0.5683464
    #> 
    #> -- sex (VI Rank: 3) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       m 0.3534041 0.2469195 0.11919360 0.5598008
    #>       f 0.2934378 0.1505113 0.04494641 0.5228192
    #> 
    #> -- copper (VI Rank: 4) --------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>    42.8 0.2626680 0.1418128 0.04426860 0.4501950
    #>    74.0 0.2798240 0.1566409 0.05285122 0.4672752
    #>     129 0.3310123 0.2146048 0.10301382 0.5322676
    #> 
    #> -- protime (VI Rank: 5) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>    10.0 0.2810872 0.1544425 0.04299654 0.4957516
    #>    10.6 0.2924100 0.1555794 0.05372951 0.5128332
    #>    11.2 0.3131481 0.1913604 0.06917838 0.5549923
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
