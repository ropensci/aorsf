
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
    #>         bili          age      ascites          ast        stage          sex 
    #>  0.017086893  0.015940821  0.006928527  0.005417795  0.004896854  0.004115441 
    #>      protime     alk.phos       hepato        edema       copper      spiders 
    #>  0.003490310  0.002760992  0.002552615  0.002449667  0.002083767  0.001198166 
    #>         trig      albumin          trt         chol 
    #>  0.000833507  0.000000000 -0.002031673 -0.002135862
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>    bili      mean        lwr      medn       upr
    #> 1:    1 0.2339613 0.01790645 0.1326649 0.8512640
    #> 2:    2 0.2848095 0.04389968 0.1814294 0.8892332
    #> 3:    3 0.3409010 0.07587147 0.2450113 0.8998751
    #> 4:    4 0.3947575 0.10778983 0.3204749 0.9189296
    #> 5:    5 0.4347783 0.14966184 0.3696747 0.9274147
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
    #>   0.80 0.2292164 0.1248095 0.05000889 0.3623311
    #>   1.40 0.2496266 0.1512690 0.06475378 0.3775441
    #>   3.52 0.3708541 0.2863213 0.18004803 0.5280146
    #> 
    #> -- age (VI Rank: 2) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   41.5 0.2697032 0.1423389 0.04636938 0.4634395
    #>   49.7 0.2945164 0.1716192 0.05345728 0.5011278
    #>   56.6 0.3259717 0.2142627 0.07571560 0.5400806
    #> 
    #> -- ascites (VI Rank: 3) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      0 0.2935145 0.1625777 0.05605702 0.5146003
    #>      1 0.4539750 0.3601905 0.26419442 0.6382479
    #> 
    #> -- ast (VI Rank: 4) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   82.5 0.2821079 0.1504670 0.04893005 0.4924931
    #>    117 0.2961175 0.1664858 0.05516819 0.5304177
    #>    153 0.3162856 0.1822462 0.06817708 0.5646494
    #> 
    #> -- stage (VI Rank: 5) --------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      1 0.2588895 0.1385436 0.04700469 0.4420654
    #>      2 0.2684688 0.1444067 0.04882304 0.4645993
    #>      3 0.2892506 0.1617974 0.05342530 0.5126863
    #>      4 0.3342546 0.2113693 0.08119181 0.5630472
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
