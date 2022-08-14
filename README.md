
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

The `orsf()` function can fit several types of ORSF ensembles. My
personal favorite is the accelerated ORSF because it has a great
combination of prediction accuracy and computational efficiency (see
[arXiv paper](https://arxiv.org/abs/2208.01129)).<sup>2</sup>

``` r
library(aorsf)

fit <- orsf(data = pbc_orsf, 
            formula = Surv(time, status) ~ . - id)
```

ORSF CONTROL VIGNETTE (not yet written).

### Inspect

Printing the output from `orsf()` will give some information and
descriptive statistics about the ensemble.

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

See
[`print.orsf_fit`](https://bcjaeger.github.io/aorsf/reference/print.orsf_fit.html)
for a description of each line in the printed output.

### Variable importance

The importance of individual variables can be estimated in three ways
using `aorsf`:

-   **negation<sup>2</sup>**: Each variable is assessed separately by
    multiplying the variable’s coefficients by -1 and then determining
    how much the model’s performance changes. The worse the model’s
    performance after negating coefficients for a given variable, the
    more important the variable.

    ``` r
    orsf_vi_negate(fit)
    #>          bili           age       ascites       protime         stage 
    #>  0.0145342780  0.0129193582  0.0059387372  0.0053657012  0.0044280058 
    #>           sex           ast         edema       spiders        hepato 
    #>  0.0031256512  0.0030735570  0.0025612975  0.0023963326  0.0015628256 
    #>          chol           trt       albumin          trig      alk.phos 
    #> -0.0001041884 -0.0009376954 -0.0018753907 -0.0020316733 -0.0022400500 
    #>      platelet 
    #> -0.0050010419
    ```

-   **permutation**: Each variable is assessed separately by randomly
    permuting the variable’s values and then determining how much the
    model’s performance changes. The worse the model’s performance after
    permuting the values of a given variable, the more important the
    variable.

    ``` r
    orsf_vi_permute(fit)
    #>          bili           age       ascites        copper         stage 
    #>  0.0109397791  0.0107314024  0.0058866431  0.0052615128  0.0041154407 
    #>       albumin        hepato       protime           ast          chol 
    #>  0.0037507814  0.0036465930  0.0022921442  0.0020316733  0.0014586372 
    #>           sex         edema       spiders      alk.phos      platelet 
    #>  0.0009376954  0.0005085385  0.0004688477  0.0000000000 -0.0007293186 
    #>           trt          trig 
    #> -0.0021358616 -0.0026047093
    ```

-   **analysis of variance (ANOVA)<sup>3</sup>**: A p-value is computed
    for each coefficient in each linear combination of variables in each
    decision tree. Importance for an individual predictor variable is
    the proportion of times a p-value for its coefficient is \< 0.01.

    ``` r
    orsf_vi_anova(fit)
    #>    ascites       bili      edema     copper    albumin        age    protime 
    #> 0.38801054 0.28122545 0.24939489 0.19932782 0.17807174 0.17757256 0.15706127 
    #>       chol      stage    spiders        ast     hepato        sex       trig 
    #> 0.14453431 0.14069149 0.13533835 0.11947845 0.11722272 0.11598746 0.10058386 
    #>   alk.phos   platelet        trt 
    #> 0.08801054 0.08376827 0.06616625
    ```

You can supply your own R function to estimate out-of-bag error when
using negation or permutation importance (see [oob
vignette](https://bcjaeger.github.io/aorsf/articles/oobag.html)).

### Partial dependence

`aorsf` can generate individual conditional expectation (ICE) and
partial dependence (PD):

-   ICE shows how the predicted value for an observation changes when a
    predictor changes. `orsf_ice()`, the function to compute ICE values
    in `aorsf`, returns ICE for all observations in the data you supply
    it. (If no data are supplied, the ORSF’s training data are used)

    ``` r
    # ICE values:
    # predicted risk for a single observation with respect to bili.

    orsf_ice_oob(fit, pd_spec = list(bili = c(1:5)))[id_row == 2]
    #>    pred_horizon id_variable id_row  bili       pred
    #>           <num>       <int>  <int> <int>      <num>
    #> 1:         1788           1      2     1 0.08949299
    #> 2:         1788           2      2     2 0.14281718
    #> 3:         1788           3      2     3 0.21558478
    #> 4:         1788           4      2     4 0.32627813
    #> 5:         1788           5      2     5 0.37660385
    ```

-   PD shows how the **expected** predicted value changes when a
    predictor changes. PD estimates the relationship between the model’s
    prediction and a predictor, adjusting for other predictors.

    ``` r
    orsf_pd_oob(fit, pd_spec = list(bili = c(1:5)))
    #>    pred_horizon  bili      mean        lwr      medn       upr
    #>           <num> <int>     <num>      <num>     <num>     <num>
    #> 1:         1788     1 0.2388058 0.01352057 0.1271245 0.8732819
    #> 2:         1788     2 0.2892796 0.03851130 0.1873985 0.9092272
    #> 3:         1788     3 0.3471000 0.07320863 0.2561256 0.9331508
    #> 4:         1788     4 0.4005475 0.10759175 0.3177516 0.9370998
    #> 5:         1788     5 0.4425184 0.14229014 0.3797550 0.9448201
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
    #>    0.80 0.2330705 0.1197286 0.04451579 0.3589746
    #>     1.4 0.2556177 0.1433992 0.06075292 0.3986335
    #>     3.5 0.3757922 0.2960614 0.16980376 0.5462387
    #> 
    #> -- age (VI Rank: 2) -----------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      42 0.2761173 0.1395697 0.04048322 0.4585127
    #>      50 0.3038888 0.1613333 0.04653317 0.5173196
    #>      57 0.3317393 0.2132545 0.06578600 0.5537216
    #> 
    #> -- ascites (VI Rank: 3) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       0 0.2977296 0.1595002 0.04905852 0.5373169
    #>       1 0.4731001 0.3866557 0.27719961 0.6610602
    #> 
    #> -- protime (VI Rank: 4) -------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>      10 0.2852584 0.1490673 0.04601975 0.5025813
    #>      11 0.2966223 0.1506002 0.05030526 0.5180246
    #>      11 0.3190773 0.1811543 0.07063537 0.5511301
    #> 
    #> -- stage (VI Rank: 5) ---------------------------
    #> 
    #>         |---------------- risk ----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       1 0.2631046 0.1353130 0.04844881 0.4482245
    #>       2 0.2726887 0.1399101 0.04656849 0.4720537
    #>       3 0.2938536 0.1608046 0.05053952 0.5175462
    #>       4 0.3421659 0.2068753 0.08564739 0.5755682
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
