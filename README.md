
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
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
#>           OOB stat value: 0.84
#>            OOB stat type: Harrell's C-statistic
#> 
#> -----------------------------------------
```

How about interpreting the fit?

-   use `orsf_vi_negate()` and `orsf_vi_anova()` for variable importance

    ``` r
    orsf_vi_negate(fit)
    #>          bili           age        copper       spiders       protime 
    #>  1.536778e-02  1.219004e-02  7.084809e-03  5.626172e-03  5.469890e-03 
    #>           sex       ascites         stage         edema           ast 
    #>  5.417795e-03  4.948948e-03  2.917274e-03  2.649361e-03  1.927485e-03 
    #>        hepato      alk.phos          trig          chol      platelet 
    #> -5.209419e-05 -1.041884e-04 -6.772244e-04 -2.865180e-03 -3.281934e-03 
    #>           trt       albumin 
    #> -3.854970e-03 -4.948948e-03
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>    bili      mean        lwr      medn       upr
    #> 1:    1 0.2378573 0.01669367 0.1259874 0.8702671
    #> 2:    2 0.2868094 0.04058178 0.1848911 0.8816975
    #> 3:    3 0.3383654 0.07527398 0.2436279 0.9040289
    #> 4:    4 0.3849085 0.10160435 0.3085406 0.9116147
    #> 5:    5 0.4250282 0.13164934 0.3626074 0.9207825
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
    #> -- bili (VI Rank: 1) ---------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   0.80 0.2333341 0.1204040 0.04807244 0.3642765
    #>   1.40 0.2528205 0.1442862 0.06114393 0.3976806
    #>   3.52 0.3644238 0.2802464 0.16067492 0.5397841
    #> 
    #> -- age (VI Rank: 2) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   41.5 0.2696389 0.1413862 0.04380715 0.4765406
    #>   49.7 0.2966209 0.1664586 0.05202837 0.5107028
    #>   56.6 0.3242722 0.2053223 0.06893931 0.5453661
    #> 
    #> -- copper (VI Rank: 3) -------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   42.8 0.2601379 0.1350286 0.04909634 0.4438070
    #>   74.0 0.2783735 0.1625457 0.05740346 0.4815901
    #>    129 0.3334568 0.2261581 0.10684394 0.5359228
    #> 
    #> -- spiders (VI Rank: 4) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      0 0.2890576 0.1529351 0.04860402 0.5119895
    #>      1 0.3369367 0.2286488 0.09168433 0.5679850
    #> 
    #> -- protime (VI Rank: 5) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   10.0 0.2798049 0.1460859 0.04749563 0.4903838
    #>   10.6 0.2924995 0.1624857 0.05138638 0.5380686
    #>   11.2 0.3140488 0.1897644 0.06859875 0.5232480
    #> 
    #>  Predicted risk at time t = 1788 for top 5 predictors
    ```

## References

Byron C. Jaeger, D. Leann Long, Dustin M. Long, Mario Sims, Jeff M.
Szychowski, Yuan-I Min, Leslie A. Mcclure, George Howard, Noah Simon
(2019). Oblique Random Survival Forests. Ann. Appl. Stat. 13(3):
1847-1883. URL <https://doi.org/10.1214/19-AOAS1261> DOI:
10.1214/19-AOAS1261
