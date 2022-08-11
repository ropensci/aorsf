
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

-   over 400 times faster than `obliqueRSF` (see Jaeger, 2019).

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

fit <- orsf(data = pbc_orsf,
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
    #>          bili           sex       protime           age         stage 
    #>  1.187747e-02  6.720150e-03  6.407585e-03  6.042926e-03  5.730360e-03 
    #>        copper       spiders       ascites         edema           ast 
    #>  5.365701e-03  3.177745e-03  2.865180e-03  2.421139e-03  2.135862e-03 
    #>        hepato          chol      alk.phos           trt       albumin 
    #>  1.041884e-03  5.209419e-05 -5.209419e-04 -2.135862e-03 -3.021463e-03 
    #>      platelet 
    #> -5.834549e-03
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>     bili      mean        lwr      medn       upr
    #>    <int>     <num>      <num>     <num>     <num>
    #> 1:     1 0.2448441 0.01284079 0.1312983 0.8747557
    #> 2:     2 0.2812156 0.02878862 0.1667878 0.8927174
    #> 3:     3 0.3202419 0.04351310 0.2171733 0.9050862
    #> 4:     4 0.3584524 0.06467319 0.2707156 0.9151960
    #> 5:     5 0.3910040 0.08913003 0.3138446 0.9204218
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
    #>    0.80 0.2405862 0.1252352 0.04384391 0.3843450
    #>     1.4 0.2568429 0.1461177 0.05419999 0.4047258
    #>     3.5 0.3412388 0.2453132 0.12975781 0.5353945
    #> 
    #> -- sex (VI Rank: 2) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       m 0.3507520 0.2415970 0.11829094 0.5454135
    #>       f 0.2908204 0.1512969 0.04511101 0.5074316
    #> 
    #> -- protime (VI Rank: 3) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      10 0.2783060 0.1510614 0.04424733 0.4966815
    #>      11 0.2900296 0.1579231 0.05008208 0.5068144
    #>      11 0.3124264 0.1818876 0.06789407 0.5116768
    #> 
    #> -- age (VI Rank: 4) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      42 0.2675298 0.1358725 0.04121726 0.4503862
    #>      50 0.2916598 0.1720650 0.04691813 0.5062555
    #>      57 0.3236244 0.2226635 0.06622152 0.5318751
    #> 
    #> -- stage (VI Rank: 5) ---------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       1 0.2513571 0.1402757 0.04384579 0.4245918
    #>       2 0.2587388 0.1458657 0.04156162 0.4391262
    #>       3 0.2829240 0.1528120 0.05249676 0.5021513
    #>       4 0.3423412 0.2281354 0.08819499 0.5622738
    #> 
    #>  Predicted risk at time t = 1788 for top 5 predictors
    ```

## Comparison to existing software

Jaeger (2022) describes `aorsf` in detail, with a summary of the
procedures used in the tree fitting algorithm and a general benchmark
comparing `aorsf` with `obliqueRSF` (and several other learners) in
terms of prediction accuracy and computational efficiency.

## References

Jaeger BC, Long DL, Long DM, Sims M, Szychowski JM, Min YI, Mcclure LA,
Howard G, Simon N. Oblique random survival forests. The Annals of
Applied Statistics. 2019 Sep;13(3):1847-83. URL:
<https://doi.org/10.1214/19-AOAS1261> DOI: 10.1214/19-AOAS1261

Jaeger BC, Welden S, Lenoir K, Speiser JL, Segar M, Pandey A, Pajewski
NM. Accelerated and interpretable oblique random survival forests. arXiv
e-prints. 2022 Aug 3:arXiv-2208. URL: <https://arxiv.org/abs/2208.01129>

## Funding

The developers of `aorsf` receive financial support from the Center for
Biomedical Informatics, Wake Forest University School of Medicine. We
also receive support from the National Center for Advancing
Translational Sciences of the National Institutes of Health under Award
Number UL1TR001420.

The content is solely the responsibility of the authors and does not
necessarily represent the official views of the National Institutes of
Health.
