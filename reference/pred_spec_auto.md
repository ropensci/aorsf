# Automatic variable values for dependence

For partial dependence and individual conditional expectations, this
function allows a variable to be considered without having to specify
what values to set the variable at. The values used are based on
quantiles for continuous variables (10th, 25th, 50th, 75th, and 90th)
and unique categories for categorical variables.

## Usage

``` r
pred_spec_auto(...)
```

## Arguments

- ...:

  names of the variables to use. These can be in quotes or not in quotes
  (see examples).

## Value

a character vector with the names

## Details

This function should only be used in the context of `orsf_pd` or
`orsf_ice` functions.

## Examples

``` r
fit <- orsf(penguins_orsf, species ~., n_tree = 5)

orsf_pd_oob(fit, pred_spec_auto(flipper_length_mm))
#> Key: <class>
#>         class flipper_length_mm      mean         lwr       medn       upr
#>        <fctr>             <num>     <num>       <num>      <num>     <num>
#>  1:    Adelie               185 0.7135339 0.052631579 0.85051447 1.0000000
#>  2:    Adelie               190 0.6938150 0.026315789 0.88479682 1.0000000
#>  3:    Adelie               197 0.6739011 0.021271930 0.90015448 1.0000000
#>  4:    Adelie               213 0.6622782 0.026315789 0.96250000 1.0000000
#>  5:    Adelie               221 0.5113873 0.017543860 0.49056604 1.0000000
#>  6: Chinstrap               185 0.2945626 0.006686137 0.12788462 1.0000000
#>  7: Chinstrap               190 0.3350056 0.006517857 0.20833333 1.0000000
#>  8: Chinstrap               197 0.3663098 0.004672897 0.18468468 1.0000000
#>  9: Chinstrap               213 0.3138748 0.007142857 0.07554359 0.9978070
#> 10: Chinstrap               221 0.2952474 0.010102960 0.18947368 0.9678195
#> 11:    Gentoo               185 0.2524412 0.016217949 0.05000000 1.0000000
#> 12:    Gentoo               190 0.3070930 0.015961538 0.07692308 1.0000000
#> 13:    Gentoo               197 0.3466848 0.004716981 0.30000000 1.0000000
#> 14:    Gentoo               213 0.5820889 0.006289308 0.65753912 1.0000000
#> 15:    Gentoo               221 0.7052264 0.009433962 0.91180654 1.0000000
```
