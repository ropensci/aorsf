# Update Forest Parameters

Update Forest Parameters

## Usage

``` r
orsf_update(object, ..., modify_in_place = FALSE, no_fit = NULL)
```

## Arguments

- object:

  (*ObliqueForest*) an oblique random forest object (see
  [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md)).

- ...:

  arguments to plug into
  [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md) that will
  be used to define the update. These arguments include:

  - `data`

  - `formula`

  - `control`

  - `weights`

  - `n_tree`

  - `n_split`

  - `n_retry`

  - `n_thread`

  - `mtry`

  - `sample_with_replacement`

  - `sample_fraction`

  - `leaf_min_events`

  - `leaf_min_obs`

  - `split_rule`

  - `split_min_events`

  - `split_min_obs`

  - `split_min_stat`

  - `pred_type`

  - `oobag_pred_horizon`

  - `oobag_eval_every`

  - `oobag_fun`

  - `importance`

  - `importance_max_pvalue`

  - `group_factors`

  - `tree_seeds`

  - `na_action`

  - `verbose_progress`

  Note that you can update `control`, but you cannot change the type of
  forest. For example, you can't go from classification to regression
  with `orsf_update`.

- modify_in_place:

  (*logical*) if `TRUE`, `object` will be modified by the inputs
  specified in `...`. Be cautious, as modification in place will
  overwrite existing data. If `FALSE` (the default), `object` will be
  copied and then the modifications will be applied to the copy, leaving
  the original `object` unmodified.

- no_fit:

  (*logical*) if `TRUE`, model fitting steps are defined and saved, but
  training is not initiated. The object returned can be directly
  submitted to
  [`orsf_train()`](https://bcjaeger.github.io/aorsf/reference/orsf.md)
  so long as `attach_data` is `TRUE`.

## Value

an `ObliqueForest` object.

## Details

There are several dynamic inputs in `orsf` with default values of
`NULL`. Specifically, these inputs are `control`, `weights`, `mtry`,
`split_rule`, `split_min_stat`, `pred_type`, `pred_horizon`,
`oobag_eval_function`, `tree_seeds`, and `oobag_eval_every`. If no
explicit value is given for these inputs in the call, they *will be
re-formed*. For example, if an initial forest includes 17 predictors,
the default `mtry` is the smallest integer that is greater than or equal
to the square root of 17, i.e., 5. Then, if you make an updated forest
with 1 less predictor and you do not explicitly say `mtry = 5`, then
`mtry` will be re-initialized in the update based on the available 16
predictors, and the resulting value of `mtry` will be 4. This is done to
avoid many potential errors that would occur if the dynamic outputs were
not re-initialized.

## Examples

``` r
if (FALSE) { # \dontrun{
# initial fit has mtry of 5
fit <- orsf(pbc_orsf, time + status ~ . -id)

# note that mtry is now 4 (see details)
fit_new <- orsf_update(fit, formula = . ~ . - edema, n_tree = 100)

# prevent dynamic updates by specifying inputs you want to freeze.
fit_newer <- orsf_update(fit_new, mtry = 2)
} # }

```
