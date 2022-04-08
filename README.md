
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aorsf

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/bcjaeger/aorsf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcjaeger/aorsf?branch=master)
[![R-CMD-check](https://github.com/bcjaeger/aorsf/workflows/R-CMD-check/badge.svg)](https://github.com/bcjaeger/aorsf/actions)
<!-- badges: end -->

The goal of `aorsf` is to fit, interpret, and make predictions with
oblique random survival forests (ORSFs). The ‘a’ in the title stands for
accelerated. So why do ORSFs need to be accelerated? Oblique decision
trees are often more accurate but slower to fit compared to their
axis-based counterparts. The issue of higher computation time is
compounded for survival decision trees, which usually require more
computing than classification or regression trees. To make ORSF more
accessible and able to engage with larger datasets, `aorsf` applies
strategies to cut down computing time without sacrificing prediction
accuracy.

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
    #>          bili           age        copper       protime       spiders 
    #>  0.0150552198  0.0102104605  0.0071369035  0.0061992082  0.0052615128 
    #>       ascites         stage           ast           sex         edema 
    #>  0.0037507814  0.0028130861  0.0026047093  0.0022400500  0.0020899691 
    #>        hepato       albumin          trig           trt      alk.phos 
    #>  0.0005209419 -0.0016149198 -0.0017712023 -0.0025005209 -0.0026568035 
    #>          chol      platelet 
    #> -0.0039070640 -0.0058345489
    ```

-   use `orsf_pd_ice()` or `orsf_pd_summary()` for individual or
    aggregated partial dependence values.

    ``` r
    orsf_pd_summary(fit, pd_spec = list(bili = c(1:5)))
    #>    bili      mean        lwr      medn       upr
    #> 1:    1 0.2376476 0.01635153 0.1256940 0.8708692
    #> 2:    2 0.2835588 0.04178780 0.1755214 0.8949604
    #> 3:    3 0.3334560 0.07231476 0.2427077 0.9133485
    #> 4:    4 0.3828746 0.10197692 0.2978942 0.9196410
    #> 5:    5 0.4246769 0.13559889 0.3548921 0.9288099
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
    #>   0.80 0.2335217 0.1199713 0.05114585 0.3770848
    #>   1.40 0.2513566 0.1450413 0.06171788 0.3921046
    #>   3.52 0.3599267 0.2732172 0.15551597 0.5235124
    #> 
    #> -- age (VI Rank: 2) ----------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   41.5 0.2696310 0.1368879 0.04219276 0.4531681
    #>   49.7 0.2977741 0.1648595 0.05141696 0.5128280
    #>   56.6 0.3290685 0.2174152 0.07230906 0.5623700
    #> 
    #> -- copper (VI Rank: 3) -------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   42.8 0.2622928 0.1338668 0.04827596 0.4526263
    #>   74.0 0.2766374 0.1574435 0.05479530 0.4745227
    #>    129 0.3308975 0.2200028 0.10231560 0.5417730
    #> 
    #> -- protime (VI Rank: 4) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>   10.0 0.2782086 0.1447316 0.04814394 0.4867040
    #>   10.6 0.2903538 0.1580250 0.05512409 0.5146609
    #>   11.2 0.3128801 0.1860102 0.06449054 0.5435381
    #> 
    #> -- spiders (VI Rank: 5) ------------------------
    #> 
    #>        |---------------- risk ----------------|
    #>  Value      Mean    Median     25th %    75th %
    #>      0 0.2889144 0.1564381 0.04787306 0.5173158
    #>      1 0.3337453 0.2158072 0.08894740 0.5493756
    #> 
    #>  Predicted risk at time t = 1788 for top 5 predictors
    ```

## References

Byron C. Jaeger, D. Leann Long, Dustin M. Long, Mario Sims, Jeff M.
Szychowski, Yuan-I Min, Leslie A. Mcclure, George Howard, Noah Simon
(2019). Oblique Random Survival Forests. Ann. Appl. Stat. 13(3):
1847-1883. URL <https://doi.org/10.1214/19-AOAS1261> DOI:
10.1214/19-AOAS1261
