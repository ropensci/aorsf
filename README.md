
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
    #>          bili           age           sex        copper       ascites 
    #>  0.0160971036  0.0120858512  0.0084392582  0.0068243384  0.0055740779 
    #>       spiders         stage         edema        hepato           ast 
    #>  0.0052094186  0.0043238175  0.0022673374  0.0021358616  0.0012502605 
    #>      platelet       albumin          chol          trig       protime 
    #>  0.0005209419  0.0003646593  0.0001562826 -0.0001562826 -0.0010939779 
    #>      alk.phos           trt 
    #> -0.0015107314 -0.0040112523
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>     bili      mean        lwr      medn       upr
    #>    <int>     <num>      <num>     <num>     <num>
    #> 1:     1 0.2309318 0.01708277 0.1168329 0.8542309
    #> 2:     2 0.2801009 0.04071749 0.1739319 0.8854798
    #> 3:     3 0.3378653 0.07326663 0.2480462 0.9002184
    #> 4:     4 0.3923115 0.10844580 0.3137873 0.9213191
    #> 5:     5 0.4341055 0.13426231 0.3711324 0.9299818
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
    #>    0.80 0.2256092 0.1153241 0.04558336 0.3588942
    #>    1.40 0.2461640 0.1313467 0.06323910 0.3871898
    #>    3.52 0.3671496 0.2847901 0.16566117 0.5423875
    #> 
    #> -- age (VI Rank: 2) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>    41.5 0.2701663 0.1369651 0.04236446 0.4578859
    #>    49.7 0.2966234 0.1673933 0.05118979 0.5077654
    #>    56.6 0.3277697 0.2096862 0.07272597 0.5657241
    #> 
    #> -- sex (VI Rank: 3) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       m 0.3485034 0.2425707 0.10799689 0.5822585
    #>       f 0.2914514 0.1429063 0.05327824 0.5329945
    #> 
    #> -- copper (VI Rank: 4) --------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>    42.8 0.2625758 0.1362109 0.04630979 0.4659081
    #>    74.0 0.2767632 0.1496957 0.05403858 0.4821049
    #>     129 0.3274266 0.2137786 0.09872470 0.5488018
    #> 
    #> -- ascites (VI Rank: 5) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       0 0.2916558 0.1506221 0.05413326 0.5223010
    #>       1 0.4545639 0.3665732 0.25809619 0.6536201
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

The developers of `aorsf` receive financial support from the Center for
Biomedical Informatics, Wake Forest University School of Medicine. We
also receive support from the National Center for Advancing
Translational Sciences of the National Institutes of Health under Award
Number UL1TR001420.

The content is solely the responsibility of the authors and does not
necessarily represent the official views of the National Institutes of
Health.
