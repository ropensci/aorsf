
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
    #>          bili           age        copper       protime       ascites 
    #>  0.0148468431  0.0131798291  0.0079183163  0.0067722442  0.0058866431 
    #>           sex         stage       spiders       albumin          chol 
    #>  0.0050010419  0.0044801000  0.0041154407  0.0034382163  0.0026568035 
    #>           ast        hepato         edema      alk.phos          trig 
    #>  0.0021879558  0.0019274849  0.0017662410  0.0017191081  0.0016149198 
    #>      platelet           trt 
    #> -0.0007293186 -0.0010418837
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>    bili      mean        lwr      medn       upr
    #> 1:    1 0.2412125 0.01670208 0.1340954 0.8821256
    #> 2:    2 0.2859373 0.04183888 0.1732189 0.9063661
    #> 3:    3 0.3350771 0.06416006 0.2333547 0.9166708
    #> 4:    4 0.3836380 0.08782332 0.3048166 0.9346609
    #> 5:    5 0.4251849 0.12446714 0.3713077 0.9424735
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
    #>   0.80 0.2368506 0.1303583 0.05151933 0.3656603
    #>   1.40 0.2546279 0.1445572 0.05960105 0.3878093
    #>   3.52 0.3613510 0.2680013 0.15962937 0.5226071
    #> 
    #> -- age (VI Rank: 2) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   41.5 0.2742174 0.1499905 0.04796971 0.4690580
    #>   49.7 0.2981615 0.1717596 0.05554041 0.5122948
    #>   56.6 0.3276850 0.2138449 0.07428326 0.5592573
    #> 
    #> -- copper (VI Rank: 3) -------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   42.8 0.2618107 0.1419154 0.05177957 0.4390637
    #>   74.0 0.2789127 0.1563250 0.05925726 0.4699412
    #>    129 0.3356129 0.2250119 0.11273457 0.5377224
    #> 
    #> -- protime (VI Rank: 4) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   10.0 0.2802231 0.1500655 0.05066870 0.4811713
    #>   10.6 0.2936040 0.1579534 0.05801067 0.5116226
    #>   11.2 0.3163518 0.1960582 0.07485445 0.5489558
    #> 
    #> -- ascites (VI Rank: 5) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      0 0.2942087 0.1592785 0.05348378 0.5189215
    #>      1 0.4672023 0.3856285 0.27823611 0.6417693
    #> 
    #>  Predicted risk at time t = 1788 for top 5 predictors
    ```

## References

Byron C. Jaeger, D. Leann Long, Dustin M. Long, Mario Sims, Jeff M.
Szychowski, Yuan-I Min, Leslie A. Mcclure, George Howard, Noah Simon
(2019). Oblique Random Survival Forests. Ann. Appl. Stat. 13(3):
1847-1883. URL <https://doi.org/10.1214/19-AOAS1261> DOI:
10.1214/19-AOAS1261
