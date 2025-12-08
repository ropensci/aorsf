# Estimate training time

Estimate training time

## Usage

``` r
orsf_time_to_train(object, n_tree_subset = NULL)
```

## Arguments

- object:

  an untrained `aorsf` object

- n_tree_subset:

  (*integer*) how many trees should be fit in order to estimate the time
  needed to train `object`. The default value is 10% of the trees
  specified in `object`. I.e., if `object` has `n_tree` of 500, then the
  default value `n_tree_subset` is 50.

## Value

a [difftime](https://rdrr.io/r/base/difftime.html) object.

## Examples

``` r
# specify but do not train the model by setting no_fit = TRUE.
object <- orsf(pbc_orsf, Surv(time, status) ~ . - id,
               n_tree = 10, no_fit = TRUE)

# approximate the time it will take to grow 10 trees
time_estimated <- orsf_time_to_train(object, n_tree_subset=1)

print(time_estimated)
#> Time difference of 0.04367828 secs

# let's see how close the approximation was
time_true_start <- Sys.time()
orsf_train(object)
time_true_stop <- Sys.time()

time_true <- time_true_stop - time_true_start

print(time_true)
#> Time difference of 0.03956723 secs

# error
abs(time_true - time_estimated)
#> Time difference of 0.004111052 secs
```
