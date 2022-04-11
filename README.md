
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
    #>           age          bili           sex        copper       protime 
    #>  0.0155761617  0.0154198791  0.0102625547  0.0090643884  0.0072410919 
    #>         stage       ascites       albumin       spiders         edema 
    #>  0.0061471140  0.0058866431  0.0048968535  0.0036986872  0.0031380546 
    #>          trig           ast        hepato          chol           trt 
    #>  0.0026047093  0.0025526151  0.0019274849  0.0001562826  0.0001041884 
    #>      alk.phos      platelet 
    #> -0.0015107314 -0.0025526151
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>    bili      mean        lwr      medn       upr
    #> 1:    1 0.2327783 0.01822227 0.1228340 0.8757152
    #> 2:    2 0.2798627 0.03933816 0.1740419 0.8969007
    #> 3:    3 0.3364899 0.06910041 0.2437527 0.9254284
    #> 4:    4 0.3929890 0.09598256 0.3157147 0.9272015
    #> 5:    5 0.4430681 0.15567318 0.3788360 0.9377603
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
    #> -- age (VI Rank: 1) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   41.5 0.2693172 0.1390511 0.04549474 0.4509366
    #>   49.7 0.2958146 0.1718434 0.05202033 0.5070320
    #>   56.6 0.3279529 0.2230427 0.07150100 0.5562976
    #> 
    #> -- bili (VI Rank: 2) ---------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   0.80 0.2286038 0.1221944 0.04762782 0.3428194
    #>   1.40 0.2475567 0.1393398 0.05955256 0.3720580
    #>   3.52 0.3647460 0.2829966 0.16185760 0.5336567
    #> 
    #> -- sex (VI Rank: 3) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      m 0.3440961 0.2421729 0.10610779 0.5679941
    #>      f 0.2926350 0.1556989 0.05196257 0.5262054
    #> 
    #> -- copper (VI Rank: 4) -------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   42.8 0.2621710 0.1437662 0.04987859 0.4348029
    #>   74.0 0.2789576 0.1609078 0.05886084 0.4653229
    #>    129 0.3311097 0.2256010 0.10596380 0.5317802
    #> 
    #> -- protime (VI Rank: 5) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   10.0 0.2791015 0.1571898 0.04683369 0.4818272
    #>   10.6 0.2913967 0.1608154 0.05454041 0.5067178
    #>   11.2 0.3134820 0.1913437 0.06815870 0.5360752
    #> 
    #>  Predicted risk at time t = 1788 for top 5 predictors
    ```

## References

Byron C. Jaeger, D. Leann Long, Dustin M. Long, Mario Sims, Jeff M.
Szychowski, Yuan-I Min, Leslie A. Mcclure, George Howard, Noah Simon
(2019). Oblique Random Survival Forests. Ann. Appl. Stat. 13(3):
1847-1883. URL <https://doi.org/10.1214/19-AOAS1261> DOI:
10.1214/19-AOAS1261
