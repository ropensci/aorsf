
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
#>  Average leaves per tree: 20
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
    #>          bili           age       protime        copper       albumin 
    #>  0.0130756408  0.0114086268  0.0098458012  0.0073973745  0.0056261721 
    #>       spiders       ascites         stage           sex           ast 
    #>  0.0050010419  0.0048447593  0.0028651802  0.0026568035  0.0019795791 
    #>         edema        hepato          trig           trt          chol 
    #>  0.0016533702  0.0015628256 -0.0007814128 -0.0015628256 -0.0026568035 
    #>      platelet 
    #> -0.0051052303
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>    bili      mean        lwr      medn       upr
    #> 1:    1 0.2361736 0.01839454 0.1180382 0.8590437
    #> 2:    2 0.2806029 0.03994244 0.1749138 0.8879428
    #> 3:    3 0.3320496 0.06504773 0.2340415 0.9022515
    #> 4:    4 0.3897615 0.11408469 0.3126960 0.9197574
    #> 5:    5 0.4354408 0.16286191 0.3630698 0.9288308
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
    #>   0.80 0.2319720 0.1160500 0.04808352 0.3702140
    #>   1.40 0.2499153 0.1359647 0.06258755 0.3931161
    #>   3.52 0.3638905 0.2809810 0.16098160 0.5329988
    #> 
    #> -- age (VI Rank: 2) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   41.5 0.2726001 0.1387599 0.04506039 0.4643828
    #>   49.7 0.2978844 0.1751892 0.05095761 0.4964155
    #>   56.6 0.3280455 0.2212188 0.06674200 0.5565002
    #> 
    #> -- protime (VI Rank: 3) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   10.0 0.2818471 0.1555631 0.04844540 0.4992442
    #>   10.6 0.2933758 0.1668372 0.05174262 0.5228040
    #>   11.2 0.3147079 0.1913758 0.06748312 0.5403061
    #> 
    #> -- copper (VI Rank: 4) -------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   42.8 0.2620289 0.1431331 0.04808758 0.4417854
    #>   74.0 0.2790941 0.1607429 0.05319370 0.4768399
    #>    129 0.3312052 0.2227699 0.09596998 0.5421644
    #> 
    #> -- albumin (VI Rank: 5) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   3.31 0.3149738 0.1839701 0.05776364 0.5530028
    #>   3.54 0.2922352 0.1571160 0.04749076 0.5216508
    #>   3.77 0.2770757 0.1478748 0.04574065 0.4912170
    #> 
    #>  Predicted risk at time t = 1788 for top 5 predictors
    ```

## References

Byron C. Jaeger, D. Leann Long, Dustin M. Long, Mario Sims, Jeff M.
Szychowski, Yuan-I Min, Leslie A. Mcclure, George Howard, Noah Simon
(2019). Oblique Random Survival Forests. Ann. Appl. Stat. 13(3):
1847-1883. URL <https://doi.org/10.1214/19-AOAS1261> DOI:
10.1214/19-AOAS1261
