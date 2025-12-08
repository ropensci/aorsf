# Inspect Forest Parameters

Printing an ORSF model tells you:

- Linear combinations: How were these identified?

- N observations: Number of rows in training data

- N events: Number of events in training data

- N trees: Number of trees in the forest

- N predictors total: Total number of columns in the predictor matrix

- N predictors per node: Number of variables used in linear combinations

- Average leaves per tree: A proxy for the depth of your trees

- Min observations in leaf: See `leaf_min_obs` in
  [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md)

- Min events in leaf: See `leaf_min_events` in
  [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md)

- OOB stat value: Out-of-bag error after fitting all trees

- OOB stat type: How was out-of-bag error computed?

- Variable importance: How was variable importance computed?

## Usage

``` r
# S3 method for class 'ObliqueForest'
print(x, ...)
```

## Arguments

- x:

  (*ObliqueForest*) an oblique random survival forest (ORSF; see
  [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md)).

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

`x`, invisibly.

## Examples

``` r
object <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 5)

print(object)
#> ---------- Oblique random survival forest
#> 
#>      Linear combinations: Accelerated Cox regression
#>           N observations: 276
#>                 N events: 111
#>                  N trees: 5
#>       N predictors total: 17
#>    N predictors per node: 5
#>  Average leaves per tree: 21.2
#> Min observations in leaf: 5
#>       Min events in leaf: 1
#>           OOB stat value: 0.76
#>            OOB stat type: Harrell's C-index
#>      Variable importance: anova
#> 
#> -----------------------------------------
```
