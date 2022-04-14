
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
[![R-CMD-check](https://github.com/bcjaeger/aorsf/workflows/R-CMD-check/badge.svg)](https://github.com/bcjaeger/aorsf/actions)
[![pkgcheck](https://github.com/bcjaeger/aorsf/workflows/pkgcheck/badge.svg)](https://github.com/%3Corg%3E/%3Crepo%3E/actions?query=workflow%3Apkgcheck)
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
    #>          bili           age       spiders        copper       ascites 
    #>  0.0145863722  0.0098978954  0.0054698896  0.0053657012  0.0041675349 
    #>         edema       protime           sex        hepato          trig 
    #>  0.0020862481  0.0016149198  0.0015628256  0.0008856012  0.0008335070 
    #>           ast      alk.phos       albumin           trt          chol 
    #>  0.0003646593 -0.0017712023 -0.0023963326 -0.0026047093 -0.0034382163 
    #>      platelet 
    #> -0.0043759116
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>    bili      mean        lwr      medn       upr
    #> 1:    1 0.2382753 0.01773748 0.1307939 0.8577912
    #> 2:    2 0.2853975 0.04275394 0.1810123 0.8854531
    #> 3:    3 0.3365038 0.07167211 0.2416212 0.9025165
    #> 4:    4 0.3866597 0.10783117 0.3041577 0.9133883
    #> 5:    5 0.4320569 0.14703665 0.3754485 0.9237173
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
    #>   0.80 0.2336123 0.1243892 0.04891279 0.3728749
    #>   1.40 0.2529711 0.1472202 0.06224278 0.4043204
    #>   3.52 0.3620186 0.2753187 0.15692651 0.5408907
    #> 
    #> -- age (VI Rank: 2) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   41.5 0.2725735 0.1390167 0.04693890 0.4488588
    #>   49.7 0.2981903 0.1660931 0.05541691 0.4978383
    #>   56.6 0.3277967 0.2120451 0.07091227 0.5494069
    #> 
    #> -- spiders (VI Rank: 3) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      0 0.2895505 0.1520606 0.04922667 0.4961561
    #>      1 0.3336321 0.2229523 0.08399896 0.5631630
    #> 
    #> -- copper (VI Rank: 4) -------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   42.8 0.2662774 0.1470582 0.04937255 0.4772178
    #>   74.0 0.2800090 0.1646179 0.05653250 0.4909421
    #>    129 0.3287680 0.2176701 0.10020571 0.5388215
    #> 
    #> -- ascites (VI Rank: 5) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      0 0.2944291 0.1537523 0.05330642 0.5204440
    #>      1 0.4536219 0.3645048 0.25408417 0.6381112
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
