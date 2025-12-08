# Variable Interactions

Use the variable interaction score described in Greenwell et al (2018).
As this method can be computationally demanding, using `n_thread=0` can
substantially reduce time needed to compute scores.

## Usage

``` r
orsf_vint(
  object,
  predictors = NULL,
  n_thread = NULL,
  verbose_progress = NULL,
  sep = ".."
)
```

## Arguments

- object:

  (*ObliqueForest*) a trained oblique random forest object (see
  [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md))

- predictors:

  (*character*) a vector of length 2 or more with names of predictors
  used by `object`. All pairwise interactions between the predictors
  will be scored. If `NULL` (the default), all predictors are used.

- n_thread:

  (*integer*) number of threads to use while growing trees, computing
  predictions, and computing importance. Default is 0, which allows a
  suitable number of threads to be used based on availability.

- verbose_progress:

  (*logical*) if `TRUE`, progress messages are printed in the console.
  If `FALSE` (the default), nothing is printed.

- sep:

  (*character*) how to separate the names of two predictors. The default
  value of `".."` returns names as `name1..name2`

## Value

a data.table with variable interaction scores and partial dependence
values.

## Details

The number of possible interactions grows exponentially based on the
number of predictors. Some caution is warranted when using large
predictor sets and it is recommended that you supply a specific vector
of predictor names to assess rather than a global search. A good
strategy is to use `n_tree = 5` to search all predictors, then pick the
top 10 interactions, get the unique predictors from them, and re-run on
just those predictors with more trees.

## References

1.  Greenwell, M B, Boehmke, C B, McCarthy, J A (2018). "A simple and
    effective model-based variable importance measure." *arXiv preprint
    arXiv:1805.04755*.

## Examples

``` r
set.seed(329)

data <- data.frame(
 x1 = rnorm(500),
 x2 = rnorm(500),
 x3 = rnorm(500)
)

data$y = with(data, expr = x1 + x2 + x3 + 1/2*x1 * x2 + x2 * x3 + rnorm(500))

forest <- orsf(data, y ~ ., n_tree = 5)

orsf_vint(forest)
#>    interaction     score          pd_values
#>         <char>     <num>             <list>
#> 1:      x2..x3 0.8021932 <data.table[25x9]>
#> 2:      x1..x2 0.5095065 <data.table[25x9]>
#> 3:      x1..x3 0.1133252 <data.table[25x9]>
```
