
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf <a href="https://bcjaeger.github.io/aorsf"><img src="man/figures/logo.png" align="right" height="138" /></a>

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
#>  Average leaves per tree: 24
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
    #>          bili           age       protime       ascites       albumin 
    #>  0.0129714524  0.0125026047  0.0085955407  0.0063554907  0.0061992082 
    #>       spiders           sex        copper         edema           ast 
    #>  0.0057303605  0.0049489477  0.0044280058  0.0027783566  0.0018753907 
    #>        hepato          trig         stage      alk.phos      platelet 
    #>  0.0018232965  0.0016670140  0.0005209419 -0.0006772244 -0.0008856012 
    #>          chol           trt 
    #> -0.0016670140 -0.0025526151
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>    bili      mean        lwr      medn       upr
    #> 1:    1 0.2376652 0.01286480 0.1331328 0.8669733
    #> 2:    2 0.2904995 0.04154464 0.1889294 0.8943053
    #> 3:    3 0.3441193 0.06285072 0.2538335 0.9146453
    #> 4:    4 0.3973732 0.10091026 0.3247072 0.9271855
    #> 5:    5 0.4424865 0.14005297 0.3807676 0.9315099
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
    #>   0.80 0.2337275 0.1267164 0.04190426 0.3716894
    #>   1.40 0.2526511 0.1452163 0.05476829 0.3998048
    #>   3.52 0.3730196 0.2841226 0.16033444 0.5673070
    #> 
    #> -- age (VI Rank: 2) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   41.5 0.2736650 0.1332843 0.04270990 0.4575449
    #>   49.7 0.2993351 0.1656678 0.04910598 0.5203980
    #>   56.6 0.3314766 0.2113831 0.06956163 0.5472662
    #> 
    #> -- protime (VI Rank: 3) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   10.0 0.2815397 0.1430703 0.04704506 0.4987604
    #>   10.6 0.2947893 0.1526181 0.05252053 0.5393637
    #>   11.2 0.3176024 0.1831842 0.06925762 0.5566000
    #> 
    #> -- ascites (VI Rank: 4) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      0 0.2954662 0.1462719 0.04778412 0.5409697
    #>      1 0.4597940 0.3742798 0.25999819 0.6551916
    #> 
    #> -- albumin (VI Rank: 5) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   3.31 0.3176012 0.1789175 0.05525473 0.5785138
    #>   3.54 0.2934144 0.1527704 0.04236004 0.5225606
    #>   3.77 0.2791547 0.1412617 0.04199867 0.4937002
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
