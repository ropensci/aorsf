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
#>  1:    Adelie               185 0.7282071 0.052631579 0.85168269 1.0000000
#>  2:    Adelie               190 0.7017577 0.026315789 0.86153274 1.0000000
#>  3:    Adelie               197 0.6804181 0.021271930 0.86337335 1.0000000
#>  4:    Adelie               213 0.6632070 0.026315789 0.96250000 1.0000000
#>  5:    Adelie               221 0.5685831 0.026315789 0.52631579 1.0000000
#>  6: Chinstrap               185 0.2958030 0.010119048 0.16203704 1.0000000
#>  7: Chinstrap               190 0.3412330 0.008363095 0.22269536 1.0000000
#>  8: Chinstrap               197 0.3782161 0.007142857 0.19369369 1.0000000
#>  9: Chinstrap               213 0.3208287 0.007142857 0.08761024 1.0000000
#> 10: Chinstrap               221 0.3198872 0.010489060 0.21052632 0.9758772
#> 11:    Gentoo               185 0.2524412 0.016217949 0.05000000 1.0000000
#> 12:    Gentoo               190 0.3057868 0.015993590 0.07692308 1.0000000
#> 13:    Gentoo               197 0.3572390 0.004716981 0.27491745 1.0000000
#> 14:    Gentoo               213 0.5727859 0.006289308 0.60000000 1.0000000
#> 15:    Gentoo               221 0.6644321 0.007783019 0.86824324 1.0000000
```
