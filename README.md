
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
[![DOI](https://zenodo.org/badge/394311897.svg)](https://zenodo.org/badge/latestdoi/394311897)
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
  hasn’t been studied as much as permutation importance. See [Jaeger,
  2022](https://arxiv.org/abs/2208.01129) for more details on this
  technique.

  ``` r

  orsf_vi_negate(fit)
  #>          bili           sex        copper           ast           age 
  #>  0.1190578208  0.0619364315  0.0290605798  0.0260108174  0.0251162396 
  #>         stage       protime         edema       ascites        hepato 
  #>  0.0237810058  0.0158443269  0.0117270641  0.0105685230  0.0092028195 
  #>       albumin          chol           trt      alk.phos       spiders 
  #>  0.0082647861  0.0041510636  0.0036548364  0.0010239241 -0.0003298163 
  #>          trig      platelet 
  #> -0.0011111508 -0.0045314656
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
  #>          bili        copper           ast           age           sex 
  #>  0.0514084384  0.0170611427  0.0142227933  0.0140274813  0.0131527430 
  #>         stage       protime       ascites         edema       albumin 
  #>  0.0119752045  0.0102865556  0.0098067817  0.0081730899  0.0080568255 
  #>        hepato          chol      alk.phos          trig       spiders 
  #>  0.0069734562  0.0032811220  0.0015862128  0.0014909643  0.0007811902 
  #>           trt      platelet 
  #> -0.0007067631 -0.0022135241
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
  #>    ascites       bili      edema        sex     copper        age        ast 
  #> 0.39107612 0.36316990 0.36316238 0.24720893 0.20547180 0.19213732 0.19029233 
  #>    albumin      stage     hepato       trig       chol    protime   alk.phos 
  #> 0.17219680 0.17068758 0.16126761 0.13379872 0.12964021 0.12659698 0.12352611 
  #>    spiders   platelet        trt 
  #> 0.11728395 0.08997135 0.07305095
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
#> -- bili (VI Rank: 1) ----------------------------
#> 
#>        |----------------- risk -----------------|
#>  Value      Mean     Median     25th %    75th %
#>   0.70 0.2074286 0.09039332 0.03827337 0.3146957
#>    1.3 0.2261739 0.10784929 0.04915971 0.3425934
#>    3.2 0.3071951 0.21242141 0.11889617 0.4358309
#> 
#> -- sex (VI Rank: 2) -----------------------------
#> 
#>        |----------------- risk -----------------|
#>  Value      Mean    Median     25th %    75th %
#>      m 0.3648659 0.2572239 0.15554270 0.5735661
#>      f 0.2479179 0.1021787 0.04161796 0.3591612
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
[arXiv paper](https://arxiv.org/abs/2208.01129). The paper

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
