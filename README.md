
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
[![R-CMD-check](https://github.com/bcjaeger/aorsf/workflows/R-CMD-check/badge.svg)](https://github.com/bcjaeger/aorsf/actions)
[![pkgcheck](https://github.com/bcjaeger/aorsf/workflows/pkgcheck/badge.svg)](https://github.com/bcjaeger/aorsf/actions?query=workflow%3Apkgcheck)
<!-- badges: end -->

`aorsf` provides optimized software to fit, interpret, and make
predictions with oblique random survival forests (ORSFs).

## Why aorsf?

-   over 500 times faster than `obliqueRSF`.

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

fit <- orsf(data_train = pbc_orsf,
            formula = Surv(time, status) ~ . - id)

print(fit)
#> ---------- Oblique random survival forest
#> 
#>           N observations: 276
#>                 N events: 111
#>                  N trees: 500
#>       N predictors total: 17
#>    N predictors per node: 5
#>  Average leaves per tree: 19
#> Min observations in leaf: 5
#>       Min events in leaf: 1
#>           OOB stat value: 0.84
#>            OOB stat type: Harrell's C-statistic
#> 
#> -----------------------------------------
```

How about interpreting the fit?

-   use `orsf_vi_negate()` and `orsf_vi_anova()` for variable importance

    ``` r
    orsf_vi_negate(fit)
    #>          bili           age       spiders       ascites        copper 
    #>  0.0156282559  0.0127109815  0.0065638675  0.0059387372  0.0046363826 
    #>       protime         edema         stage        hepato          trig 
    #>  0.0045842884  0.0035386837  0.0022921442  0.0013023547  0.0010418837 
    #>           sex           ast           trt      alk.phos      platelet 
    #>  0.0007293186  0.0006251302 -0.0011981663 -0.0018753907 -0.0019795791 
    #>          chol       albumin 
    #> -0.0023963326 -0.0047405710
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>    bili      mean        lwr      medn       upr
    #> 1:    1 0.2352732 0.01833523 0.1196051 0.8675898
    #> 2:    2 0.2821136 0.04112874 0.1714480 0.8919186
    #> 3:    3 0.3371132 0.06911118 0.2545069 0.9177142
    #> 4:    4 0.3922269 0.10551068 0.3141783 0.9270324
    #> 5:    5 0.4376877 0.13935792 0.3658866 0.9337528
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
    #> -- bili (VI Rank: 1) ---------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   0.80 0.2302554 0.1169780 0.04856717 0.3521751
    #>   1.40 0.2486047 0.1329087 0.06193197 0.3859956
    #>   3.52 0.3669522 0.2779187 0.16387346 0.5419349
    #> 
    #> -- age (VI Rank: 2) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   41.5 0.2720260 0.1353055 0.04497034 0.4635398
    #>   49.7 0.2961534 0.1684676 0.05106990 0.5080294
    #>   56.6 0.3249653 0.2054728 0.06901834 0.5609686
    #> 
    #> -- spiders (VI Rank: 3) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      0 0.2900883 0.1547913 0.04929566 0.4958708
    #>      1 0.3327155 0.2041412 0.08623750 0.5486246
    #> 
    #> -- ascites (VI Rank: 4) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      0 0.2942298 0.1605606 0.05226003 0.5278228
    #>      1 0.4616882 0.3752462 0.25852339 0.6415898
    #> 
    #> -- copper (VI Rank: 5) -------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   42.8 0.2634105 0.1427394 0.04614312 0.4569345
    #>   74.0 0.2789849 0.1600635 0.05747877 0.4679921
    #>    129 0.3305907 0.2252363 0.09389927 0.5383923
    #> 
    #>  Predicted risk at time t = 1788 for top 5 predictors
    ```

## References

Byron C. Jaeger, D. Leann Long, Dustin M. Long, Mario Sims, Jeff M.
Szychowski, Yuan-I Min, Leslie A. Mcclure, George Howard, Noah Simon
(2019). Oblique Random Survival Forests. Ann. Appl. Stat. 13(3):
1847-1883. URL <https://doi.org/10.1214/19-AOAS1261> DOI:
10.1214/19-AOAS1261

## Funding

This software receives financial support from the Center for Biomedical
Informatics, Wake Forest School of Medicine.
