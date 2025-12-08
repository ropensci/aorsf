# Introduction to aorsf

This article covers core features of the `aorsf` package.

``` r
library(aorsf)
```

## Background

The oblique random forest (RF) is an extension of the traditional
(axis-based) RF. Instead of using a single variable to split data and
grow new branches, trees in the oblique RF use a weighted combination of
multiple variables.

## Oblique RFs for survival, classification, and regression

The purpose of `aorsf` (‘a’ is short for accelerated) is to provide a
unifying framework to fit oblique RFs that can scale adequately to large
data sets. The fastest algorithms available in the package are used by
default because they often have equivalent prediction accuracy to more
computational approaches.

The center piece of `aorsf` is the
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) function.
In the initial versions of `aorsf`, the
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) function
only fit `o`blique `r`andom `s`urvival `f`orests, but now it allows for
classification, regression, and survival forests. (I may introduce an
`orf()` function in the future if the name
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) is
misleading to users.)

For classification, we fit an oblique RF to predict penguin species
using `penguin` data from the magnificent `palmerpenguins` [R
package](https://allisonhorst.github.io/palmerpenguins/)

``` r
# An oblique classification RF
penguin_fit <- orsf(data = penguins_orsf, formula = species ~ .)

penguin_fit
#> ---------- Oblique random classification forest
#> 
#>      Linear combinations: Accelerated Logistic regression
#>           N observations: 333
#>                N classes: 3
#>                  N trees: 500
#>       N predictors total: 7
#>    N predictors per node: 3
#>  Average leaves per tree: 5.522
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
bill_fit <- orsf(data = penguins_orsf, formula = bill_length_mm ~ .)

bill_fit
#> ---------- Oblique random regression forest
#> 
#>      Linear combinations: Accelerated Linear regression
#>           N observations: 333
#>                  N trees: 500
#>       N predictors total: 7
#>    N predictors per node: 3
#>  Average leaves per tree: 49.806
#> Min observations in leaf: 5
#>           OOB stat value: 0.81
#>            OOB stat type: RSQ
#>      Variable importance: anova
#> 
#> -----------------------------------------
```

My personal favorite is the oblique survival RF with accelerated Cox
regression because it has a great combination of prediction accuracy and
computational efficiency (see [JCGS
paper](https://doi.org/10.1080/10618600.2023.2231048)). Here, we predict
mortality risk following diagnosis of primary biliary cirrhosis:

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
#>  Average leaves per tree: 21.6
#> Min observations in leaf: 5
#>       Min events in leaf: 1
#>           OOB stat value: 0.77
#>            OOB stat type: Harrell's C-index
#>      Variable importance: anova
#> 
#> -----------------------------------------
```

you may notice that the first input of `aorsf` is `data`. This is a
design choice that makes it easier to use `orsf` with pipes (i.e., `%>%`
or `|>`). For instance,

``` r

library(dplyr)

pbc_fit <- pbc_orsf |> 
 select(-id) |> 
 orsf(formula = Surv(time, status) ~ .,
      n_tree = 5)
```

## Interpretation

`aorsf` includes several functions dedicated to interpretation of ORSFs,
both through estimation of partial dependence and variable importance.

### Variable importance

There are multiple methods to compute variable importance, and each can
be applied to any type of oblique forest.

- To compute *negation* importance, ORSF multiplies each coefficient of
  that variable by -1 and then re-computes the out-of-sample (sometimes
  referred to as out-of-bag) accuracy of the ORSF model.

  ``` r

  orsf_vi_negate(pbc_fit)
  #>          bili           age        copper           ast           sex 
  #>  0.1468851774  0.0606952129  0.0246435580  0.0224269123  0.0175587328 
  #>          trig      alk.phos       protime         edema          chol 
  #>  0.0096895007  0.0093198869  0.0086039712  0.0006382134 -0.0015687436 
  #>       ascites      platelet        hepato       spiders           trt 
  #> -0.0060269468 -0.0102280228 -0.0108549805 -0.0113883544 -0.0201827916 
  #>         stage       albumin 
  #> -0.0221462608 -0.0224072750
  ```

- You can also compute variable importance using *permutation*, a more
  classical approach that noises up a predictor and then assigned the
  resulting degradation in prediction accuracy to be the importance of
  that predictor.

  ``` r

  orsf_vi_permute(penguin_fit)
  #>    bill_length_mm flipper_length_mm     bill_depth_mm            island 
  #>      0.1718853563      0.1019742036      0.0740207434      0.0707193201 
  #>       body_mass_g               sex              year 
  #>      0.0636020709      0.0179647016      0.0007314716
  ```

- A faster alternative to permutation and negation importance is ANOVA
  importance, which computes the proportion of times each variable
  obtains a low p-value (p \< 0.01) while the forest is grown.

  ``` r

  orsf_vi_anova(bill_fit)
  #>           species               sex            island flipper_length_mm 
  #>        0.34766885        0.20967742        0.11680487        0.08913484 
  #>       body_mass_g     bill_depth_mm              year 
  #>        0.07659396        0.05951482        0.01587443
  ```

### Partial dependence (PD)

Partial dependence (PD) shows the expected prediction from a model as a
function of a single predictor or multiple predictors. The expectation
is marginalized over the values of all other predictors, giving
something like a multivariable adjusted estimate of the model’s
prediction.

For more on PD, see the
[vignette](https://docs.ropensci.org/aorsf/articles/pd.html)

### Individual conditional expectations (ICE)

Unlike partial dependence, which shows the expected prediction as a
function of one or multiple predictors, individual conditional
expectations (ICE) show the prediction for an individual observation as
a function of a predictor.

For more on ICE, see the
[vignette](https://docs.ropensci.org/aorsf/articles/pd.html#individual-conditional-expectations-ice)

## What about the original ORSF?

The original ORSF (i.e., `obliqueRSF`) used `glmnet` to find linear
combinations of inputs. `aorsf` allows users to implement this approach
using the `orsf_control_survival(method = 'net')` function:

``` r

orsf_net <- orsf(data = pbc_orsf, 
                 formula = Surv(time, status) ~ . - id, 
                 control = orsf_control_survival(method = 'net'))
```

`net` forests fit a lot faster than the original ORSF function in
`obliqueRSF`. However, `net` forests are still much slower than `cph`
ones.

## aorsf and other machine learning software

The unique feature of `aorsf` is its fast algorithms to fit ORSF
ensembles. `RLT` and `obliqueRSF` both fit oblique random survival
forests, but `aorsf` does so faster. `ranger` and `randomForestSRC` fit
survival forests, but neither package supports oblique splitting.
`obliqueRF` fits oblique random forests for classification and
regression, but not survival. `PPforest` fits oblique random forests for
classification but not survival.

Note: The default prediction behavior for `aorsf` models is to produce
predicted risk at a specific prediction horizon, which is not the
default for `ranger` or `randomForestSRC`. I think this will change in
the future, as computing time independent predictions with `aorsf` could
be helpful.

## Learning more

`aorsf` began as a dedicated package for oblique random survival
forests, and so most papers published so far have focused on survival
analysis and risk prediction. However, the routines for regression and
classification oblique RFs in `aorsf` have high overlap with the
survival ones.

- See [orsf](https://docs.ropensci.org/aorsf/reference/orsf.html) for
  more details on oblique random survival forests.

- see the [JCGS](https://doi.org/10.1080/10618600.2023.2231048) paper
  for more details on algorithms used specifically by `aorsf`.
