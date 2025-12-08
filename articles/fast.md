# Tips to speed up computation

``` r
library(aorsf)
```

## Go faster

Analyses can slow to a crawl when models need hours to run. In this
article you will find a few tricks to prevent this bottleneck when using
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md).

## Don’t specify a `control`

The default `control` for
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) is `NULL`
because, if unspecified,
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) will pick
the fastest possible `control` for you depending on the type of forest
being grown. The default `control` run-time compared to other approaches
can be striking. For example:

``` r

time_fast <- system.time(
 expr = orsf(pbc_orsf, 
             formula = time+status~. -id, 
             n_tree = 5)
)

time_net <- system.time(
 expr = orsf(pbc_orsf, 
             formula = time+status~. -id, 
             control = orsf_control_survival(method = 'net'), 
             n_tree = 5)
)

# unspecified control is much faster
time_net['elapsed'] / time_fast['elapsed']
#> elapsed 
#>    47.9
```

## Use `n_thread`

The `n_thread` argument uses multi-threading to run `aorsf` functions in
parallel when possible. If you know how many threads you want, e.g. you
want exactly 5, set `n_thread = 5`. If you aren’t sure how many threads
you have available but want to use a feasible amount, using
`n_thread = 0` (the default) tells `aorsf` to do that for you.

``` r

# automatically pick number of threads based on amount available

orsf(pbc_orsf, 
     formula = time+status~. -id, 
     n_tree = 5,
     n_thread = 0)
```

Note: sometimes multi-threading is not possible. For example, because R
is a single threaded language, multi-threading cannot be applied when
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) needs to
call R functions from C++, which occurs when a customized R function is
used to find linear combination of variables or compute prediction
accuracy.

## Do less

There are some inputs in
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) that can
be adjusted to make it run faster:

- set `n_retry` to `0`

- set `oobag_pred_type` to `'none'`

- set `importance` to `'none'`

- increase `split_min_events`, `split_min_obs`, `leaf_min_events`, or
  `leaf_min_obs` to make trees stop growing sooner

- increase `split_min_stat` to enforce more strict requirements for
  growing deeper trees.

Applying these tips:

``` r

orsf(pbc_orsf, 
     formula = time+status~., 
     n_thread = 0, 
     n_tree = 5, 
     n_retry = 0,
     oobag_pred_type = 'none', 
     importance = 'none',
     split_min_events = 20, 
     leaf_min_events = 10,
     split_min_stat = 10)
```

While modifying these inputs can make
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) run
faster, they can also impact prediction accuracy.

## Show progress

Setting `verbose_progress = TRUE` doesn’t make anything run faster, but
it can help make it *feel* like things are running less slow.

``` r

verbose_fit <- orsf(pbc_orsf, 
                    formula = time+status~. -id, 
                    n_tree = 5, 
                    verbose_progress = TRUE)
#> Growing trees: 100%. 
#> Computing predictions: 100%.
```

## Don’t wait. Estimate!

Instead of running a model and hoping it will be fast, you can estimate
how long a specification of that model will take by using
`no_fit = TRUE` in the call to
[`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md).

``` r

fit_spec <- orsf(pbc_orsf, 
                 formula = time+status~. -id, 
                 control = orsf_control_survival(method = 'net'), 
                 n_tree = 2000,
                 no_fit = TRUE)

# how much time it takes to estimate training time:
system.time(
 time_est <- orsf_time_to_train(fit_spec, n_tree_subset = 5)
)
#>    user  system elapsed 
#>   0.294   0.000   0.294

# the estimated training time:
time_est
#> Time difference of 117.1861 secs
```
