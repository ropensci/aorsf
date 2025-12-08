# Scale input data

These functions are exported so that users may access internal routines
that are used to scale inputs when
[orsf_control_cph](https://bcjaeger.github.io/aorsf/reference/orsf_control_cph.md)
is used.

## Usage

``` r
orsf_scale_cph(x_mat, w_vec = NULL)

orsf_unscale_cph(x_mat)
```

## Arguments

- x_mat:

  (*numeric matrix*) a matrix with values to be scaled or unscaled. Note
  that `orsf_unscale_cph` will only accept `x_mat` inputs that have an
  attribute containing transform values, which are added automatically
  by `orsf_scale_cph`.

- w_vec:

  (*numeric vector*) an optional vector of weights. If no weights are
  supplied (the default), all observations will be equally weighted. If
  supplied, `w_vec` must have length equal to `nrow(x_mat)`.

## Value

the scaled or unscaled `x_mat`.

## Details

The data are transformed by first subtracting the mean and then
multiplying by the scale. An inverse transform can be completed using
`orsf_unscale_cph` or by dividing each column by the corresponding scale
and then adding the mean.

The values of means and scales are stored in an attribute of the output
returned by `orsf_scale_cph` (see examples)

## Examples

``` r
x_mat <- as.matrix(pbc_orsf[, c('bili', 'age', 'protime')])

head(x_mat)
#>   bili      age protime
#> 1 14.5 58.76523    12.2
#> 2  1.1 56.44627    10.6
#> 3  1.4 70.07255    12.0
#> 4  1.8 54.74059    10.3
#> 5  3.4 38.10541    10.9
#> 7  1.0 55.53457     9.7

x_scaled <- orsf_scale_cph(x_mat)

head(x_scaled)
#>             bili        age    protime
#> [1,]  3.77308887  1.0412574  1.9694656
#> [2,] -0.75476469  0.7719344 -0.1822316
#> [3,] -0.65339483  2.3544852  1.7005035
#> [4,] -0.51823502  0.5738373 -0.5856748
#> [5,]  0.02240421 -1.3581657  0.2212116
#> [6,] -0.78855464  0.6660494 -1.3925613

attributes(x_scaled) # note the transforms attribute
#> $dim
#> [1] 276   3
#> 
#> $dimnames
#> $dimnames[[1]]
#> NULL
#> 
#> $dimnames[[2]]
#> [1] "bili"    "age"     "protime"
#> 
#> 
#> $transforms
#>           mean     scale
#> [1,]  3.333696 0.3378995
#> [2,] 49.799661 0.1161396
#> [3,] 10.735507 1.3448108
#> 

x_unscaled <- orsf_unscale_cph(x_scaled)

head(x_unscaled)
#>      bili      age protime
#> [1,] 14.5 58.76523    12.2
#> [2,]  1.1 56.44627    10.6
#> [3,]  1.4 70.07255    12.0
#> [4,]  1.8 54.74059    10.3
#> [5,]  3.4 38.10541    10.9
#> [6,]  1.0 55.53457     9.7

# numeric difference in x_mat and x_unscaled should be practically 0
max(abs(x_mat - x_unscaled))
#> [1] 8.881784e-16
```
