
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf <a href="https://docs.ropensci.org/aorsf/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
[![R-CMD-check](https://github.com/bcjaeger/aorsf/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/aorsf/actions/)
[![Status at rOpenSci Software Peer
Review](https://badges.ropensci.org/532_status.svg)](https://github.com/ropensci/software-review/issues/532/)
<a href="https://joss.theoj.org/papers/10.21105/joss.04705"><img src="https://joss.theoj.org/papers/10.21105/joss.04705/status.svg"></a>
[![CRAN
status](https://www.r-pkg.org/badges/version/aorsf)](https://CRAN.R-project.org/package=aorsf)
[![DOI](https://zenodo.org/badge/394311897.svg)](https://zenodo.org/doi/10.5281/zenodo.7116854)
<!-- badges: end -->

Fit, interpret, and make predictions with oblique random survival
forests (ORSFs).

## Why aorsf?

- Hundreds of times faster than `obliqueRSF`.<sup>1</sup>

- Accurate predictions for censored outcomes.<sup>2</sup>

- Negation importance, a novel technique to estimate variable importance
  for ORSFs.<sup>2</sup>

- Intuitive API with formula based interface.

- Extensive input checks and informative error messages.

## Installation

You can install `aorsf` from CRAN using

``` r
install.packages("aorsf")
```

You can install the development version of aorsf from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/aorsf")
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
[JCGS
paper](https://doi.org/10.1080/10618600.2023.2231048)).<sup>2</sup>

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
#>  Average leaves per tree: 10
#> Min observations in leaf: 5
#>       Min events in leaf: 1
#>           OOB stat value: 0.83
#>            OOB stat type: Harrell's C-statistic
#>      Variable importance: anova
#> 
#> -----------------------------------------
```

- See
  [print.orsf_fit](https://docs.ropensci.org/aorsf/reference/print.orsf_fit.html)
  for a description of each line in the printed output.

- See [orsf
  examples](https://docs.ropensci.org/aorsf/reference/orsf.html#examples)
  for more details on controlling ORSF ensemble fits and using them in
  prediction modeling workflows.

### Variable importance

The importance of individual variables can be estimated in three ways
using `aorsf`:

- **negation**<sup>2</sup>: Each variable is assessed separately by
  multiplying the variable’s coefficients by -1 and then determining how
  much the model’s performance changes. The worse the model’s
  performance after negating coefficients for a given variable, the more
  important the variable. This technique is promising b/c it does not
  require permutation and it emphasizes variables with larger
  coefficients in linear combinations, but it is also relatively new and
  hasn’t been studied as much as permutation importance. See
  Jaeger, (2023) for more details on this technique.

  ``` r

  orsf_vi_negate(fit)
  #>          bili           sex        copper         stage           age 
  #>  0.1162463868  0.0517905362  0.0375565841  0.0240450064  0.0239056901 
  #>           ast       protime        hepato         edema       ascites 
  #>  0.0191083400  0.0158014897  0.0139536512  0.0119264604  0.0100865906 
  #>       albumin          chol       spiders           trt          trig 
  #>  0.0085394443  0.0037903802  0.0030727468  0.0020617896  0.0018361632 
  #>      alk.phos      platelet 
  #>  0.0006586211 -0.0044967624
  ```

- **permutation**: Each variable is assessed separately by randomly
  permuting the variable’s values and then determining how much the
  model’s performance changes. The worse the model’s performance after
  permuting the values of a given variable, the more important the
  variable. This technique is flexible, intuitive, and frequently used.
  It also has several [known
  limitations](https://christophm.github.io/interpretable-ml-book/feature-importance.html#disadvantages-9)

  ``` r

  orsf_vi_permute(fit)
  #>          bili        copper           age         stage           sex 
  #>  0.0523994364  0.0187964038  0.0152246586  0.0115192591  0.0110127557 
  #>           ast        hepato         edema       ascites       albumin 
  #>  0.0100104477  0.0082889176  0.0079183046  0.0077834483  0.0070642325 
  #>       protime          trig          chol       spiders      alk.phos 
  #>  0.0066513097  0.0015656325  0.0014474560  0.0006015308  0.0001369292 
  #>           trt      platelet 
  #> -0.0013984860 -0.0022427356
  ```

- **analysis of variance (ANOVA)**<sup>3</sup>: A p-value is computed
  for each coefficient in each linear combination of variables in each
  decision tree. Importance for an individual predictor variable is the
  proportion of times a p-value for its coefficient is \< 0.01. This
  technique is very efficient computationally, but may not be as
  effective as permutation or negation in terms of selecting signal over
  noise variables. See [Menze,
  2011](https://link.springer.com/chapter/10.1007/978-3-642-23783-6_29)
  for more details on this technique.

  ``` r

  orsf_vi_anova(fit)
  #>       bili    ascites      edema     copper      stage        sex        age 
  #> 0.48778004 0.44943820 0.41677872 0.31865585 0.26675095 0.26458616 0.25448430 
  #>        ast     hepato    albumin       chol    protime       trig    spiders 
  #> 0.21743929 0.19945726 0.18191604 0.15240328 0.15076561 0.13709677 0.11833550 
  #>   alk.phos   platelet        trt 
  #> 0.10113636 0.06302021 0.05019305
  ```

You can supply your own R function to estimate out-of-bag error when
using negation or permutation importance (see [oob
vignette](https://docs.ropensci.org/aorsf/articles/oobag.html))

### Partial dependence (PD)

Partial dependence (PD) shows the expected prediction from a model as a
function of a single predictor or multiple predictors. The expectation
is marginalized over the values of all other predictors, giving
something like a multivariable adjusted estimate of the model’s
prediction.

The summary function, `orsf_summarize_uni()`, computes PD for as many
variables as you ask it to, using sensible values.

``` r

orsf_summarize_uni(fit, n_variables = 2)
#> 
#> -- bili (VI Rank: 1) ---------------------------
#> 
#>        |---------------- Risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>   0.70 0.1986719 0.1044026 0.04354701 0.2968290
#>    1.3 0.2132847 0.1210276 0.05245387 0.3208855
#>    3.2 0.2883814 0.1917119 0.11951296 0.4147258
#> 
#> -- sex (VI Rank: 2) ----------------------------
#> 
#>        |---------------- Risk ----------------|
#>  Value      Mean    Median     25th %    75th %
#>      m 0.3394141 0.2313787 0.13762339 0.5311308
#>      f 0.2390067 0.1112093 0.04782891 0.3773551
#> 
#>  Predicted risk at time t = 1826.25 for top 2 predictors
```

For more on PD, see the
[vignette](https://docs.ropensci.org/aorsf/articles/pd.html)

### Individual conditional expectations (ICE)

Unlike partial dependence, which shows the expected prediction as a
function of one or multiple predictors, individual conditional
expectations (ICE) show the prediction for an individual observation as
a function of a predictor.

For more on ICE, see the
[vignette](https://docs.ropensci.org/aorsf/articles/pd.html#individual-conditional-expectations-ice)

## Comparison to existing software

Comparisons between `aorsf` and existing software are presented in our
[JCGS paper](https://doi.org/10.1080/10618600.2023.2231048). The paper:

- describes `aorsf` in detail with a summary of the procedures used in
  the tree fitting algorithm

- runs a general benchmark comparing `aorsf` with `obliqueRSF` and
  several other learners

- reports prediction accuracy and computational efficiency of all
  learners.

- runs a simulation study comparing variable importance techniques with
  ORSFs, axis based RSFs, and boosted trees.

- reports the probability that each variable importance technique will
  rank a relevant variable with higher importance than an irrelevant
  variable.

A more hands-on comparison of `aorsf` and other R packages is provided
in [orsf
examples](https://docs.ropensci.org/aorsf/reference/orsf.html#tidymodels)

## References

1.  Jaeger BC, Long DL, Long DM, Sims M, Szychowski JM, Min YI, Mcclure
    LA, Howard G, Simon N. Oblique random survival forests. *Annals of
    applied statistics* 2019 Sep; 13(3):1847-83. DOI:
    10.1214/19-AOAS1261

2.  Jaeger BC, Welden S, Lenoir K, Speiser JL, Segar MW, Pandey A,
    Pajewski NM. Accelerated and interpretable oblique random survival
    forests. *Journal of Computational and Graphical Statistics*
    Published online 08 Aug 2023. DOI: 10.1080/10618600.2023.2231048

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
