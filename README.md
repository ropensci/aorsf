
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
    #>          bili           age        copper           sex       ascites 
    #>  0.0143259012  0.0135965826  0.0082829756  0.0078662221  0.0050010419 
    #>         stage       spiders       protime           ast         edema 
    #>  0.0041675349  0.0040112523  0.0036986872  0.0036465930  0.0027833180 
    #>          chol      alk.phos        hepato       albumin          trig 
    #>  0.0007293186  0.0006251302  0.0005730360  0.0004688477 -0.0017712023 
    #>      platelet           trt 
    #> -0.0032298395 -0.0036986872
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>     bili      mean        lwr      medn       upr
    #>    <int>     <num>      <num>     <num>     <num>
    #> 1:     1 0.2469777 0.01345288 0.1315400 0.8857421
    #> 2:     2 0.2839531 0.02911864 0.1684011 0.8986682
    #> 3:     3 0.3241288 0.04961583 0.2261766 0.9148388
    #> 4:     4 0.3645008 0.07139964 0.2834042 0.9188526
    #> 5:     5 0.3974148 0.09856039 0.3182883 0.9282819
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
    #>    0.80 0.2426995 0.1308204 0.05059714 0.3869034
    #>     1.4 0.2594904 0.1421002 0.06048081 0.4252482
    #>     3.5 0.3457362 0.2538523 0.13447505 0.5257414
    #> 
    #> -- age (VI Rank: 2) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      42 0.2659744 0.1379398 0.04009427 0.4327782
    #>      50 0.2950065 0.1664274 0.05262147 0.4928079
    #>      57 0.3286207 0.2179738 0.07502315 0.5466582
    #> 
    #> -- copper (VI Rank: 3) --------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      43 0.2630432 0.1471669 0.04848619 0.4517837
    #>      74 0.2799984 0.1593967 0.06055610 0.4846709
    #>     129 0.3257924 0.2135642 0.09899796 0.5250845
    #> 
    #> -- sex (VI Rank: 4) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       m 0.3538255 0.2450266 0.12591346 0.5722112
    #>       f 0.2907292 0.1571660 0.05143278 0.5015086
    #> 
    #> -- ascites (VI Rank: 5) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       0 0.2924870 0.1591568 0.05274239 0.5167522
    #>       1 0.4295291 0.3452233 0.21235610 0.6302404
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

Jaeger BC, Welden S, Lenoir K, Speiser JL, Segar MW, Pandey A, Pajewski
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
