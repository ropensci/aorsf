# Print ORSF summary

Print ORSF summary

## Usage

``` r
# S3 method for class 'orsf_summary_uni'
print(x, n_variables = NULL, ...)
```

## Arguments

- x:

  an object of class 'orsf_summary'

- n_variables:

  The number of variables to print

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

invisibly, `x`

## Examples

``` r
object <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 25)

smry <- orsf_summarize_uni(object, n_variables = 2)

print(smry)
#> 
#> -- ascites (VI Rank: 1) -------------------------
#> 
#>         |---------------- Risk ----------------|
#>   Value      Mean    Median     25th %    75th %
#>  <char>     <num>     <num>      <num>     <num>
#>       0 0.3232406 0.2049484 0.03946059 0.5767793
#>       1 0.4742986 0.4411374 0.23183711 0.6975218
#> 
#> -- bili (VI Rank: 2) ----------------------------
#> 
#>         |---------------- Risk ----------------|
#>   Value      Mean    Median     25th %    75th %
#>  <char>     <num>     <num>      <num>     <num>
#>    0.60 0.2512724 0.1460509 0.02993850 0.3969644
#>    0.80 0.2565434 0.1485614 0.03020138 0.3999826
#>    1.40 0.2765232 0.1754551 0.05600655 0.4348171
#>    3.52 0.4018822 0.3283903 0.18810057 0.5898275
#>    7.25 0.4987987 0.4546319 0.32170399 0.6635148
#> 
#>  Predicted risk at time t = 1788 for top 2 predictors 
```
