
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

-   Hundreds of times faster than `obliqueRSF` (see Jaeger, 2019).

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
#>                ORSF type: Accelerated
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
#>      Variable importance: anova
#> 
#> -----------------------------------------
```

How about interpreting the fit?

-   use `orsf_vi_negate()` and `orsf_vi_anova()` for variable importance

    ``` r
    orsf_vi_negate(fit)
    #>          bili           age       ascites       protime       albumin 
    #>  0.0171910815  0.0149510315  0.0057303605  0.0052615128  0.0034903105 
    #>         stage       spiders      platelet        copper         edema 
    #>  0.0034382163  0.0033340279  0.0030214628  0.0027609919  0.0022425307 
    #>           ast        hepato          trig           trt      alk.phos 
    #>  0.0021358616  0.0006251302 -0.0002604709 -0.0005730360 -0.0013544488 
    #>           sex          chol 
    #> -0.0017191081 -0.0027088977
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>     bili      mean        lwr      medn       upr
    #>    <int>     <num>      <num>     <num>     <num>
    #> 1:     1 0.2365249 0.01101961 0.1321714 0.8670140
    #> 2:     2 0.2871905 0.03974545 0.1744483 0.8973028
    #> 3:     3 0.3427409 0.07308052 0.2560156 0.9234381
    #> 4:     4 0.3940760 0.10322698 0.3167858 0.9324083
    #> 5:     5 0.4348478 0.12844698 0.3807357 0.9380762
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
    #>    0.80 0.2307432 0.1265528 0.04429212 0.3632950
    #>     1.4 0.2528437 0.1414677 0.05762821 0.3833136
    #>     3.5 0.3714271 0.2887484 0.16208629 0.5438425
    #> 
    #> -- age (VI Rank: 2) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      42 0.2720845 0.1315823 0.04113686 0.4559970
    #>      50 0.3001501 0.1662907 0.04352686 0.5166783
    #>      57 0.3317869 0.2131719 0.06890813 0.5587767
    #> 
    #> -- ascites (VI Rank: 3) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median    25th %    75th %
    #>  <char>     <num>     <num>     <num>     <num>
    #>       0 0.2955811 0.1572229 0.0468181 0.5181648
    #>       1 0.4641456 0.3752439 0.2604392 0.6542257
    #> 
    #> -- protime (VI Rank: 4) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      10 0.2821317 0.1483479 0.04539954 0.4939123
    #>      11 0.2935816 0.1557434 0.05169954 0.5268243
    #>      11 0.3167064 0.1849317 0.06648985 0.5525649
    #> 
    #> -- albumin (VI Rank: 5) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>     3.3 0.3186478 0.1883154 0.05615637 0.5750781
    #>     3.5 0.2938590 0.1562017 0.04646367 0.5291128
    #>     3.8 0.2788236 0.1419260 0.04421349 0.4806785
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
