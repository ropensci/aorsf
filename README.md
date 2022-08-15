
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

## What is an oblique decision tree?

Decision trees are developed by splitting a set of training data into
two new subsets, with the goal of having more similarity within the new
subsets than between them. The splitting process is repeated on
resulting subsets of data until a stopping criterion is met.

When the new subsets of data are formed based on a single predictor, the
decision tree is said to be *axis-based* because the splits of the data
appear perpendicular to the axis of the predictor. When linear
combinations of variables are used instead of a single variable, the
tree is *oblique* because the splits of the data are neither parallel
nor at a right angle to the axis.

**Figure**: Decision trees for classification with axis-based splitting
(left) and oblique splitting (right). Cases are orange squares; controls
are purple circles. Both trees partition the predictor space defined by
variables X1 and X2, but the oblique splits do a better job of
separating the two classes.

<img src="man/figures/tree_axis_v_oblique.png" width="100%" />

## Examples

The `orsf()` function can fit several types of ORSF ensembles. My
personal favorite is the accelerated ORSF because it has a great
combination of prediction accuracy and computational efficiency (see
[arXiv paper](https://arxiv.org/abs/2208.01129)).<sup>2</sup>

``` r
library(aorsf)

set.seed(329730)

index_train <- sample(nrow(pbc_orsf), 150) 

pbc_orsf_train <- pbc_orsf[index_train, ]
pbc_orsf_test <- pbc_orsf[-index_train, ]

fit <- orsf(data = pbc_orsf_train, 
            formula = Surv(time, status) ~ . - id,
            oobag_pred_horizon = 365.25 * 5)
```

ORSF CONTROL VIGNETTE (not yet written).

### Inspect

Printing the output from `orsf()` will give some information and
descriptive statistics about the ensemble.

``` r
fit
#> ---------- Oblique random survival forest
#> 
#>      Linear combinations: Accelerated
#>           N observations: 150
#>                 N events: 52
#>                  N trees: 500
#>       N predictors total: 17
#>    N predictors per node: 5
#>  Average leaves per tree: 12
#> Min observations in leaf: 5
#>       Min events in leaf: 1
#>           OOB stat value: 0.83
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

-   **negation**<sup>2</sup>: Each variable is assessed separately by
    multiplying the variable’s coefficients by -1 and then determining
    how much the model’s performance changes. The worse the model’s
    performance after negating coefficients for a given variable, the
    more important the variable.

    ``` r
    orsf_vi_negate(fit)
    #>          bili           sex           age       ascites         edema 
    #>  0.0152354571  0.0138504155  0.0132568263  0.0059358924  0.0051286110 
    #>         stage      alk.phos        hepato       protime        copper 
    #>  0.0023743569  0.0011871785  0.0007914523  0.0005935892 -0.0001978631 
    #>           ast       albumin           trt      platelet          chol 
    #> -0.0005935892 -0.0021764939 -0.0041551247 -0.0043529877 -0.0051444400
    ```

-   **permutation**: Each variable is assessed separately by randomly
    permuting the variable’s values and then determining how much the
    model’s performance changes. The worse the model’s performance after
    permuting the values of a given variable, the more important the
    variable.

    ``` r
    orsf_vi_permute(fit)
    #>         bili      ascites          age          sex        edema      albumin 
    #>  0.010091017  0.007716660  0.007320934  0.006727345  0.004337159  0.003561535 
    #>        stage      protime       hepato         chol      spiders       copper 
    #>  0.003165809  0.002770083  0.002176494  0.001187178 -0.001780768 -0.002176494 
    #>     platelet          trt         trig 
    #> -0.002572220 -0.004155125 -0.004946577
    ```

-   **analysis of variance (ANOVA)**<sup>3</sup>: A p-value is computed
    for each coefficient in each linear combination of variables in each
    decision tree. Importance for an individual predictor variable is
    the proportion of times a p-value for its coefficient is \< 0.01.

    ``` r
    orsf_vi_anova(fit)
    #>    ascites       bili      edema        sex        age     copper      stage 
    #> 0.35231788 0.33216374 0.31401592 0.22045995 0.19044776 0.18155620 0.16907605 
    #>        ast     hepato    albumin       chol       trig    protime    spiders 
    #> 0.14183124 0.13736655 0.12611012 0.11461988 0.10847044 0.10697115 0.08802817 
    #>   alk.phos   platelet        trt 
    #> 0.07943094 0.06150342 0.04411765
    ```

You can supply your own R function to estimate out-of-bag error when
using negation or permutation importance (see [oob
vignette](https://bcjaeger.github.io/aorsf/articles/oobag.html)).

### Predict

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

    orsf_ice_oob(fit, pred_spec = list(bili = c(1:5)))[id_row == 2]
    #>    pred_horizon id_variable id_row  bili       pred
    #>           <num>       <int>  <int> <int>      <num>
    #> 1:      1826.25           1      2     1 0.05397791
    #> 2:      1826.25           2      2     2 0.09138356
    #> 3:      1826.25           3      2     3 0.12458702
    #> 4:      1826.25           4      2     4 0.17812030
    #> 5:      1826.25           5      2     5 0.25564427
    ```

-   PD shows how the **expected** predicted value changes when a
    predictor changes. PD estimates the relationship between the model’s
    prediction and a predictor, adjusting for other predictors.

    ``` r
    orsf_pd_oob(fit, pred_spec = list(bili = c(1:5)))
    #>    pred_horizon  bili      mean        lwr       medn       upr
    #>           <num> <int>     <num>      <num>      <num>     <num>
    #> 1:      1826.25     1 0.2075896 0.01389732 0.09063976 0.7998756
    #> 2:      1826.25     2 0.2352634 0.02628113 0.12935779 0.8152149
    #> 3:      1826.25     3 0.2750782 0.04254451 0.18877830 0.8371582
    #> 4:      1826.25     4 0.3302680 0.08806724 0.24827784 0.8441472
    #> 5:      1826.25     5 0.3846734 0.14808075 0.29926304 0.8562432
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
    #> -- bili (VI Rank: 1) -----------------------------
    #> 
    #>         |----------------- risk -----------------|
    #>   Value      Mean     Median     25th %    75th %
    #>  <char>     <num>      <num>      <num>     <num>
    #>    0.70 0.2021325 0.08822768 0.03120946 0.3037762
    #>     1.3 0.2150109 0.09563599 0.04159316 0.3138399
    #>     3.2 0.2835993 0.19877449 0.09747990 0.4360340
    #> 
    #> -- sex (VI Rank: 2) ------------------------------
    #> 
    #>         |----------------- risk -----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       m 0.3639036 0.2690788 0.15537470 0.5608909
    #>       f 0.2379669 0.1025245 0.03534165 0.3794900
    #> 
    #> -- age (VI Rank: 3) ------------------------------
    #> 
    #>         |----------------- risk -----------------|
    #>   Value      Mean     Median     25th %    75th %
    #>  <char>     <num>      <num>      <num>     <num>
    #>      41 0.2307407 0.09229274 0.03250819 0.3657784
    #>      50 0.2478440 0.10323803 0.03373931 0.4087894
    #>      57 0.2767813 0.14199710 0.04893782 0.4888879
    #> 
    #> -- ascites (VI Rank: 4) --------------------------
    #> 
    #>         |----------------- risk -----------------|
    #>   Value      Mean    Median    25th %    75th %
    #>  <char>     <num>     <num>     <num>     <num>
    #>       0 0.2504481 0.1128150 0.0368409 0.5055972
    #>       1 0.4329801 0.3504737 0.2724964 0.5838057
    #> 
    #> -- edema (VI Rank: 5) ----------------------------
    #> 
    #>         |----------------- risk -----------------|
    #>   Value      Mean    Median     25th %    75th %
    #>  <char>     <num>     <num>      <num>     <num>
    #>       0 0.2395980 0.1047356 0.03534165 0.5055972
    #>     0.5 0.3213759 0.2013491 0.10819771 0.5723427
    #>       1 0.4408370 0.3740121 0.27777420 0.6475894
    #> 
    #>  Predicted risk at time t = 1826.25 for top 5 predictors
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
