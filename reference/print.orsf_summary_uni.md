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
#>       0 0.3238634 0.2174355 0.03903864 0.5747091
#>       1 0.4746833 0.4410430 0.23299544 0.6975218
#> 
#> -- bili (VI Rank: 2) ----------------------------
#> 
#>         |---------------- Risk ----------------|
#>   Value      Mean    Median     25th %    75th %
#>  <char>     <num>     <num>      <num>     <num>
#>    0.60 0.2511524 0.1515735 0.02920573 0.3902603
#>    0.80 0.2567416 0.1558418 0.02895238 0.3989856
#>    1.40 0.2786348 0.1789260 0.06166630 0.4321895
#>    3.52 0.4052499 0.3440759 0.19964396 0.5861602
#>    7.25 0.4992038 0.4559516 0.32237422 0.6606566
#> 
#>  Predicted risk at time t = 1788 for top 2 predictors 
```
