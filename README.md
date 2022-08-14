
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf <a href="https://bcjaeger.github.io/aorsf"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
[![R-CMD-check](https://github.com/bcjaeger/aorsf/workflows/R-CMD-check/badge.svg)](https://github.com/bcjaeger/aorsf/actions)
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
#>      Linear combinations: Accelerated
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
    #>          bili           age       protime       spiders       ascites 
    #>  0.0189101896  0.0156282559  0.0073973745  0.0051573244  0.0044280058 
    #>        copper         stage           ast          trig           sex 
    #>  0.0040633465  0.0034903105  0.0025526151  0.0025005209  0.0018753907 
    #>         edema        hepato      platelet      alk.phos          chol 
    #>  0.0013978607  0.0007293186 -0.0006772244 -0.0009376954 -0.0014065430 
    #>           trt 
    #> -0.0019274849
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>     bili      mean        lwr      medn       upr
    #>    <int>     <num>      <num>     <num>     <num>
    #> 1:     1 0.2356302 0.01203499 0.1233501 0.8728360
    #> 2:     2 0.2872993 0.03491402 0.1809715 0.8968723
    #> 3:     3 0.3428164 0.05744697 0.2435341 0.9165324
    #> 4:     4 0.3911216 0.09348188 0.3115675 0.9336801
    #> 5:     5 0.4323540 0.12792062 0.3634836 0.9436788
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
    #>    0.80 0.2301513 0.1216147 0.04465412 0.3713734
    #>     1.4 0.2521101 0.1355369 0.05720828 0.4022439
    #>     3.5 0.3693030 0.2789399 0.15395771 0.5777261
    #> 
    #> -- age (VI Rank: 2) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      42 0.2698680 0.1360924 0.03826626 0.4375958
    #>      50 0.2984951 0.1646925 0.04567883 0.5237952
    #>      57 0.3317490 0.2179783 0.06907194 0.5597332
    #> 
    #> -- protime (VI Rank: 3) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      10 0.2803568 0.1514527 0.04432300 0.5159335
    #>      11 0.2936656 0.1588985 0.05149954 0.5392990
    #>      11 0.3162593 0.1932949 0.06516967 0.5436588
    #> 
    #> -- spiders (VI Rank: 4) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       0 0.2906128 0.1567162 0.04526985 0.5137760
    #>       1 0.3336593 0.2091174 0.08381937 0.5628989
    #> 
    #> -- ascites (VI Rank: 5) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       0 0.2948807 0.1581709 0.04695869 0.5352667
    #>       1 0.4622080 0.3779279 0.25267155 0.6586023
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
