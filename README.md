
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

-   Hundreds of times faster than `obliqueRSF`.<sup>1</sup>

-   Fast and accurate predictions for censored outcomes.<sup>2</sup>

-   negation importance, a novel technique to estimate variable
    importance for ORSFs.<sup>2</sup>

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

The `orsf()` function is used to fit ORSF ensembles:

``` r
library(aorsf)

fit <- orsf(data = pbc_orsf, 
            formula = Surv(time, status) ~ . - id)
```

The default routine to fit ORSF ensembles is the ‘accelerated’ ORSF - an
algorithm based on Newton Raphson scoring that does very well in
benchmarks of prediction accuracy and computational
efficiency.<sup>2</sup> In addition to the accelerated ORSF, `aorsf` can
fit a broad range of ORSF ensembles (see ORSF CONTROL VIGNETTE (not yet
written)).

### Inspect

Printing the output from `orsf()` will give some descriptive statistics
about the ensemble.

``` r
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

### Variable importance

The importance of individual variables can be estimated in three ways
using `aorsf`:

-   **negation**: Each variable is assessed separately by multiplying
    the variable’s coefficients by -1 and then determining how much the
    model’s performance changes. The worse the model’s performance after
    negating coefficients for a given variable, the more important the
    variable.

    ``` r
    orsf_vi_negate(fit)
    #>           age          bili        copper         stage       protime 
    #>  0.0151594082  0.0133361117  0.0065638675  0.0060950198  0.0060429256 
    #>       albumin           sex       spiders       ascites           ast 
    #>  0.0059387372  0.0058345489  0.0049489477  0.0048447593  0.0035424047 
    #>         edema        hepato          chol          trig           trt 
    #>  0.0019498110  0.0006772244 -0.0014065430 -0.0015107314 -0.0019795791 
    #>      platelet      alk.phos 
    #> -0.0026568035 -0.0032819337
    ```

-   **permutation**: Each variable is assessed separately by randomly
    permuting the variable’s values and then determining how much the
    model’s performance changes. The worse the model’s performance after
    permuting the values of a given variable, the more important the
    variable.

    ``` r
    orsf_vi_permute(fit)
    #>          bili           age        copper       albumin         stage 
    #>  1.604501e-02  1.302355e-02  5.261513e-03  3.959158e-03  3.907064e-03 
    #>       protime       ascites       spiders           sex          trig 
    #>  3.021463e-03  2.917274e-03  9.897895e-04  4.167535e-04  1.562826e-04 
    #>         edema        hepato      platelet          chol           trt 
    #> -4.713284e-05 -1.041884e-04 -1.562826e-04 -4.167535e-04 -6.772244e-04 
    #>           ast      alk.phos 
    #> -1.198166e-03 -1.354449e-03
    ```

-   **analysis of variance (ANOVA)<sup>3</sup>**: A p-value is computed
    for each coefficient in each linear combination of variables in each
    decision tree. Importance for an individual predictor variable is
    the proportion of times a p-value for its coefficient is \< 0.01.

    ``` r
    orsf_vi_anova(fit)
    #>    ascites       bili      edema     copper        age    albumin    protime 
    #> 0.35348226 0.28289811 0.24968033 0.18991641 0.18409387 0.16945107 0.15829608 
    #>      stage        ast       chol    spiders        sex     hepato       trig 
    #> 0.13969986 0.13060480 0.12707469 0.12549740 0.11944046 0.11162362 0.10188777 
    #>   alk.phos   platelet        trt 
    #> 0.09502618 0.07333506 0.05134680
    ```

You can also supply your own R function to estimate out-of-bag error
when using negation or permutation importance (see [oob
vignette](https://bcjaeger.github.io/aorsf/articles/oobag.html)).

### Partial dependence

`aorsf` can generate individual conditional expectation (ICE) and
partial dependence:

-   ICE is the expected predicted value of an ORSF ensemble for an
    individual observation.

-   partial dependence is a multi-variable adjusted expected predicted
    value of an ORSF ensemble.

    ``` r
    orsf_pd(fit, pd_spec = list(bili = c(1:5)))
    #>     bili      mean        lwr      medn       upr
    #>    <int>     <num>      <num>     <num>     <num>
    #> 1:     1 0.2338035 0.01359849 0.1221437 0.8637697
    #> 2:     2 0.2864137 0.03772655 0.1806327 0.8944082
    #> 3:     3 0.3434524 0.06803086 0.2455590 0.9050007
    #> 4:     4 0.3949240 0.09451781 0.3250008 0.9312628
    #> 5:     5 0.4372363 0.13161986 0.3747099 0.9377383
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
    #> -- age (VI Rank: 1) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      42 0.2711144 0.1371126 0.04265667 0.4737754
    #>      50 0.3007558 0.1661609 0.04900273 0.5097726
    #>      57 0.3319072 0.2089503 0.07061881 0.5643859
    #> 
    #> -- bili (VI Rank: 2) ----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>    0.80 0.2291821 0.1170113 0.04592342 0.3688073
    #>     1.4 0.2496643 0.1388448 0.06020376 0.4054563
    #>     3.5 0.3705226 0.2927490 0.16222596 0.5511317
    #> 
    #> -- copper (VI Rank: 3) --------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      43 0.2674542 0.1443824 0.04631652 0.4825262
    #>      74 0.2831247 0.1603985 0.05519124 0.5099932
    #>     129 0.3336233 0.2247753 0.10327836 0.5373395
    #> 
    #> -- stage (VI Rank: 4) ---------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       1 0.2613092 0.1346674 0.04656766 0.4589368
    #>       2 0.2715204 0.1358854 0.04618129 0.4782684
    #>       3 0.2942587 0.1557578 0.05312238 0.5204041
    #>       4 0.3386524 0.2051842 0.08445636 0.5724707
    #> 
    #> -- protime (VI Rank: 5) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      10 0.2833402 0.1498676 0.04555345 0.5049473
    #>      11 0.2960497 0.1602747 0.05321297 0.5378991
    #>      11 0.3155702 0.1887370 0.06624194 0.5559068
    #> 
    #>  Predicted risk at time t = 1788 for top 5 predictors
    ```

## Comparison to existing software

Jaeger (2022) describes `aorsf` in detail, with a summary of the
procedures used in the tree fitting algorithm and a general benchmark
comparing `aorsf` with `obliqueRSF` (and several other learners) in
terms of prediction accuracy and computational efficiency.

## References

1.  Jaeger BC, Long DL, Long DM, Sims M, Szychowski JM, Min YI, Mcclure
    LA, Howard G, Simon N. Oblique random survival forests. *Annals of
    applied statistics* 2019 Sep; 13(3):1847-83. DOI:
    10.1214/19-AOAS1261

2.  Jaeger BC, Welden S, Lenoir K, Speiser JL, Segar MW, Pandey A,
    Pajewski NM. Accelerated and interpretable oblique random survival
    forests. *arXiv e-prints* 2022 Aug; arXiv-2208. URL:
    <https://arxiv.org/abs/2208.01129>

3.  Menze BH, Kelm BM, Splitthoff DN, Koethe U, Hamprecht FA. On oblique
    random forests. *Joint European Conference on Machine Learning and
    Knowledge Discovery in Databases* 2011 Sep 4; pp. 453-469. DOI:
    10.1007/978-3-642-23783-6_29

## Funding

The developers of `aorsf` receive financial support from the Center for
Biomedical Informatics, Wake Forest University School of Medicine. We
also receive support from the National Center for Advancing
Translational Sciences of the National Institutes of Health under Award
Number UL1TR001420.

The content is solely the responsibility of the authors and does not
necessarily represent the official views of the National Institutes of
Health.
