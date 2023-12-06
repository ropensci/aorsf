
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

Fit, interpret, and make predictions with oblique random forests (RFs).

## Why aorsf?

- Fast and versatile tools for oblique RFs.<sup>1</sup>

- Accurate predictions.<sup>2</sup>

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

## Get started

``` r
library(aorsf)
library(tidyverse)
```

`aorsf` fits several types of oblique RFs with the `orsf()` function,
including survival, classification, and regression RFs. My personal
favorite is the oblique survival RF with accelerated Cox regression
because it has a great combination of prediction accuracy and
computational efficiency (see [JCGS
paper](https://doi.org/10.1080/10618600.2023.2231048)).<sup>2</sup>

For classification, we fit an oblique RF to predict penguin species:

``` r
# An oblique classification RF
penguin_fit <- orsf(data = penguins_orsf,
                    n_tree = 5, 
                    formula = species ~ .)

penguin_fit
#> ---------- Oblique random classification forest
#> 
#>      Linear combinations: Accelerated Logistic regression
#>           N observations: 333
#>                N classes: 3
#>                  N trees: 5
#>       N predictors total: 7
#>    N predictors per node: 3
#>  Average leaves per tree: 4.8
#> Min observations in leaf: 5
#>           OOB stat value: 1.00
#>            OOB stat type: AUC-ROC
#>      Variable importance: anova
#> 
#> -----------------------------------------
```

For regression, we use the same data but predict bill length of
penguins:

``` r
# An oblique regression RF
bill_fit <- orsf(data = penguins_orsf, 
                 n_tree = 5, 
                 formula = bill_length_mm ~ .)

bill_fit
#> ---------- Oblique random regression forest
#> 
#>      Linear combinations: Accelerated Linear regression
#>           N observations: 333
#>                  N trees: 5
#>       N predictors total: 7
#>    N predictors per node: 3
#>  Average leaves per tree: 50
#> Min observations in leaf: 5
#>           OOB stat value: 0.69
#>            OOB stat type: RSQ
#>      Variable importance: anova
#> 
#> -----------------------------------------
```

For survival, we predict mortality risk following diagnosis of primary
biliary cirrhosis:

``` r
# An oblique survival RF
pbc_fit <- orsf(data = pbc_orsf, 
                n_tree = 5,
                formula = Surv(time, status) ~ . - id)

pbc_fit
#> ---------- Oblique random survival forest
#> 
#>      Linear combinations: Accelerated Cox regression
#>           N observations: 276
#>                 N events: 111
#>                  N trees: 5
#>       N predictors total: 17
#>    N predictors per node: 5
#>  Average leaves per tree: 20.2
#> Min observations in leaf: 5
#>       Min events in leaf: 1
#>           OOB stat value: 0.75
#>            OOB stat type: Harrell's C-index
#>      Variable importance: anova
#> 
#> -----------------------------------------
```

## What does “oblique” mean?

Decision trees are grown by splitting a set of training data into
non-overlapping subsets, with the goal of having more similarity within
the new subsets than between them. When subsets are created with a
single predictor, the decision tree is *axis-based* because the subset
boundaries are perpendicular to the axis of the predictor. When linear
combinations (i.e., a weighted sum) of variables are used instead of a
single variable, the tree is *oblique* because the boundaries are
neither parallel nor perpendicular to the axis.

**Figure**: Decision trees for classification with axis-based splitting
(left) and oblique splitting (right). Cases are orange squares; controls
are purple circles. Both trees partition the predictor space defined by
variables X1 and X2, but the oblique splits do a better job of
separating the two classes.

<img src="man/figures/tree_axis_v_oblique.png" width="100%" />

So, how does this difference translate to real data, and how does it
impact random forests comprising hundreds of axis-based or oblique
trees? We will demonstrate this using the `penguin` data from the
magnificent `palmerpenguins` R package. We will also use this function
to make several plots:

``` r
plot_decision_surface <- function(predictions, title, grid){
 
 # this is not a general function for plotting
 # decision surfaces. It just helps to minimize 
 # copying and pasting of code.
 
 colnames(predictions) <- levels(penguins_orsf$species)
 
 class_preds <- bind_cols(grid, predictions) %>%
  pivot_longer(cols = c(Adelie,
                        Chinstrap,
                        Gentoo)) %>%
  group_by(flipper_length_mm, bill_length_mm) %>%
  arrange(desc(value)) %>%
  slice(1)
 
 cols <- c("darkorange", "purple", "cyan4")

 ggplot(class_preds, aes(bill_length_mm, flipper_length_mm)) +
  geom_contour_filled(aes(z = value, fill = name),
                      alpha = .25) +
  geom_point(data = penguins_orsf,
             aes(color = species, shape = species),
             size = 2,
             alpha = 0.8) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(x = "Bill length, mm",
       y = "Flipper length, mm") +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = '') + 
  labs(title = title)
 
}
```

We also use a grid of points for plotting decision surfaces:

``` r
grid <- expand_grid(

 flipper_length_mm = seq(min(penguins_orsf$flipper_length_mm),
                     max(penguins_orsf$flipper_length_mm),
                  len = 200),
 bill_length_mm = seq(min(penguins_orsf$bill_length_mm),
                      max(penguins_orsf$bill_length_mm),
                      len = 200)
)
```

We use `orsf` with `mtry=1` to fit axis-based trees:

``` r
fit_axis_tree <- penguins_orsf %>% 
 orsf(species ~ bill_length_mm + flipper_length_mm,
      n_tree = 1,
      mtry = 1,
      tree_seeds = 106760)
```

Next we use `orsf_update` to copy and modify the original model,
expanding it to fit an oblique tree by using `mtry=2` instead of
`mtry=1`, and to include 500 trees instead of 1:

``` r
fit_axis_forest <- fit_axis_tree %>% 
 orsf_update(n_tree = 500)

fit_oblique_tree <- fit_axis_tree %>% 
 orsf_update(mtry = 2)

fit_oblique_forest <- fit_oblique_tree %>% 
 orsf_update(n_tree = 500)
```

And now we have all we need to visualize decision surfaces using
predictions from these four fits:

``` r
preds <- list(fit_axis_tree,
              fit_axis_forest,
              fit_oblique_tree,
              fit_oblique_forest) %>% 
 map(predict, new_data = grid, pred_type = 'prob')

titles <- c("Axis-based tree",
            "Axis-based forest",
            "Oblique tree",
            "Oblique forest")

plots <- map2(preds, titles,  
              .f = plot_decision_surface, 
              grid = grid)
```

**Figure**: Axis-based and oblique decision surfaces from a single tree
and an ensemble of 500 trees. Axis-based trees have boundaries
perpendicular to predictor axes, whereas oblique trees can have
boundaries that are neither parallel nor perpendicular to predictor
axes. Axis-based forests tend to have ‘step-function’ decision
boundaries, while oblique forests tend to have smooth decision
boundaries.

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

## Variable importance

The importance of individual predictor variables can be estimated in
three ways using `aorsf` and can be used on any type of oblique RF.
Also, variable importance functions always return a named character
vector

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
  orsf_vi_negate(pbc_fit)
  #>         bili          sex       copper      protime         trig        stage 
  #>  0.123886962  0.079314182  0.033690082  0.024993393  0.024824543  0.022596171 
  #>          age      ascites         chol      albumin          trt          ast 
  #>  0.022230055  0.022102660  0.021148146  0.019735757  0.016349041  0.015422072 
  #>        edema     platelet      spiders       hepato     alk.phos 
  #>  0.011281582  0.009236892 -0.004852893 -0.008039998 -0.016934342
  ```

- **permutation**: Each variable is assessed separately by randomly
  permuting the variable’s values and then determining how much the
  model’s performance changes. The worse the model’s performance after
  permuting the values of a given variable, the more important the
  variable. This technique is flexible, intuitive, and frequently used.
  It also has several [known
  limitations](https://christophm.github.io/interpretable-ml-book/feature-importance.html#disadvantages-9)

  ``` r
  orsf_vi_permute(penguin_fit)
  #>    bill_length_mm            island flipper_length_mm       body_mass_g 
  #>       0.143264349       0.087622456       0.083339890       0.075864339 
  #>     bill_depth_mm               sex              year 
  #>       0.069297784       0.016607005       0.002510418
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
  orsf_vi_anova(bill_fit)
  #>           species               sex            island flipper_length_mm 
  #>        0.61764706        0.29166667        0.15302191        0.12403101 
  #>     bill_depth_mm       body_mass_g              year 
  #>        0.07142857        0.06106870        0.02564103
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
orsf_summarize_uni(pbc_fit, n_variables = 2)
#> 
#> -- bili (VI Rank: 1) -----------------------------
#> 
#>         |----------------- Risk -----------------|
#>   Value      Mean     Median     25th %    75th %
#>  <char>     <num>      <num>      <num>     <num>
#>    0.80 0.1962502 0.05782313 0.01521944 0.3064935
#>    1.40 0.2184447 0.06462585 0.01541984 0.3333333
#>    3.55 0.3613502 0.27833749 0.10859877 0.6000000
#> 
#> -- sex (VI Rank: 2) ------------------------------
#> 
#>         |----------------- Risk -----------------|
#>   Value      Mean    Median     25th %    75th %
#>  <char>     <num>     <num>      <num>     <num>
#>       m 0.3329003 0.2222222 0.01636962 0.5740741
#>       f 0.2692446 0.1055556 0.01541984 0.4000000
#> 
#>  Predicted risk at time t = 1788 for top 2 predictors
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
  oblique survival RFs, axis based survival RFs, and boosted trees.

- reports the probability that each variable importance technique will
  rank a relevant variable with higher importance than an irrelevant
  variable.

<!-- A more hands-on comparison of `aorsf` and other R packages is provided in [orsf examples](https://docs.ropensci.org/aorsf/reference/orsf.html#tidymodels) -->

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

The developers of `aorsf` received financial support from the Center for
Biomedical Informatics, Wake Forest University School of Medicine. We
also received support from the National Center for Advancing
Translational Sciences of the National Institutes of Health under Award
Number UL1TR001420.

The content is solely the responsibility of the authors and does not
necessarily represent the official views of the National Institutes of
Health.
