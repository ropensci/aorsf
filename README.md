
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf <a href="https://bcjaeger.github.io/aorsf/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
[![R-CMD-check](https://github.com/bcjaeger/aorsf/workflows/R-CMD-check/badge.svg)](https://github.com/bcjaeger/aorsf/actions/)
[![Status at rOpenSci Software Peer
Review](https://badges.ropensci.org/532_status.svg)](https://github.com/ropensci/software-review/issues/532/)
<a href="https://joss.theoj.org/papers/414871f081cd8449007d671a7f7f7c3a"><img src="https://joss.theoj.org/papers/414871f081cd8449007d671a7f7f7c3a/status.svg"></a>
<!-- badges: end -->

Fit, interpret, and make predictions with oblique random survival
forests (ORSFs).

## Why aorsf?

-   Hundreds of times faster than `obliqueRSF`.<sup>1</sup>

-   Accurate predictions for censored outcomes.<sup>2</sup>

-   Negation importance, a novel technique to estimate variable
    importance for ORSFs.<sup>2</sup>

-   Intuitive API with formula based interface.

-   Extensive input checks and informative error messages.

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

-   See
    [print.orsf_fit](https://bcjaeger.github.io/aorsf/reference/print.orsf_fit.html)
    for a description of each line in the printed output.

-   See [orsf
    examples](https://bcjaeger.github.io/aorsf/reference/orsf.html#examples)
    for more details on controlling ORSF ensemble fits and using them in
    prediction modeling workflows.

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
using negation or permutation importance. This feature is experimental
and may be changed in the future (see [oob
vignette](https://bcjaeger.github.io/aorsf/articles/oobag.html))

### Partial dependence (PD)

Partial dependence (PD) shows the expected prediction from a model as a
function of a single predictor or multiple predictors. The expectation
is marginalized over the values of all other predictors, giving
something like a multivariable adjusted estimate of the model’s
prediction.

For more on PD, see the
[vignette](https://bcjaeger.github.io/aorsf/articles/pd.html)

### Individual conditional expectations (ICE)

Unlike partial dependence, which shows the expected prediction as a
function of one or multiple predictors, individual conditional
expectations (ICE) show the prediction for an individual observation as
a function of a predictor.

For more on ICE, see the
[vignette](https://bcjaeger.github.io/aorsf/articles/pd.html#individual-conditional-expectations-ice)

## Comparison to existing software

Comparisons between `aorsf` and existing software are presented in our
[arXiv paper](https://arxiv.org/abs/2208.01129/). The paper

-   describes `aorsf` in detail with a summary of the procedures used in
    the tree fitting algorithm

-   runs a general benchmark comparing `aorsf` with `obliqueRSF` and
    several other learners

-   reports prediction accuracy and computational efficiency of all
    learners.

-   runs a simulation study comparing variable importance techniques
    with ORSFs, axis based RSFs, and boosted trees.

-   reports the probability that each variable importance technique will
    rank a relevant variable with higher importance than an irrelevant
    variable.

A more hands-on comparison of `aorsf` and other R packages is provided
in [orsf
examples](https://bcjaeger.github.io/aorsf/reference/orsf.html#tidymodels)

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
