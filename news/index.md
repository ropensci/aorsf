# Changelog

## aorsf 0.1.6

CRAN release: 2025-12-11

- Corrected documentation indicating principal component analysis was a
  built-in option and fixed an issue where `method="fast"` was not
  giving the expected error inside of `orsf_control` functions (see
  <https://github.com/ropensci/aorsf/pull/79>). Thank you
  [@emilyriederer](https://github.com/emilyriederer)!

- added `n_predictor_drop` to
  [`orsf_vs()`](https://bcjaeger.github.io/aorsf/reference/orsf_vs.md).
  Dropping one predictor at a time makes
  [`orsf_vs()`](https://bcjaeger.github.io/aorsf/reference/orsf_vs.md)
  slow for data with hundreds of predictors. Using a larger value for
  `n_predictor_drop` helps speed this up. The default value of
  `n_predictor_drop` is 1 to maintain backward compatibility.

- `orsf` no longer throws hard errors if `leaf_min_obs` or
  `leaf_min_events` exceed the highest admissible value. Instead, a
  warning is returned and the input value is replaced by the highest
  admissible value.

## aorsf 0.1.5

CRAN release: 2024-05-30

- fixed an issue where omitting NA values would cause an error in
  regression forests.

## aorsf 0.1.4

CRAN release: 2024-05-03

- `orsf_vs` now returns a column that contains non-reference coded
  variable names (see <https://github.com/ropensci/aorsf/pull/52>).

- `orsf_vs` no longer throws an error when `n_predictor_min = 1` is used
  (see <https://github.com/ropensci/aorsf/pull/58>).

- `orsf_summarize_uni` now allows specification of a class to summarize
  for oblique classification forests (see
  <https://github.com/ropensci/aorsf/pull/57>).

- fixed an issue where `orsf` would throw an uninformative error when
  all predictors were categorical (see
  <https://github.com/ropensci/aorsf/pull/56>)

- oblique random forests can now compute out-of-bag predictions on
  modified versions of their training data (see
  <https://github.com/ropensci/aorsf/pull/54>)

- Setting `oobag_pred_type` to `'none'` when growing a forest no longer
  necessitates the specification of `pred_type` when calling `predict`
  later (see <https://github.com/ropensci/aorsf/pull/48>).

- Setting `sample_fraction` to 1 will no longer result in empty
  `oobag_rows` in the forest object (this would cause R to crash when
  the forest was passed to C++; see
  <https://github.com/ropensci/aorsf/pull/48>)

- Re-worked the creation and maintenance of `oobag_denom` in C++
  routines (see <https://github.com/ropensci/aorsf/pull/48>).

- Restricted mean survival time is now used for `pred_type = 'time'`
  instead of median survival time (See
  <https://github.com/ropensci/aorsf/pull/46>).

## aorsf 0.1.3

CRAN release: 2024-01-22

- minor changes to partial dependence vignette to resolve code
  sanitization errors.

## aorsf 0.1.2

CRAN release: 2024-01-15

- Allowed option `"time"` for `pred_type` in `predict` and partial
  dependence to predict survival time (see
  <https://github.com/ropensci/aorsf/issues/37>).

- Added
  [`pred_spec_auto()`](https://bcjaeger.github.io/aorsf/reference/pred_spec_auto.md)
  for more convenient specification of variables for partial dependence.

- Partial dependence now runs much faster with multiple threads.

- Added
  [`orsf_vint()`](https://bcjaeger.github.io/aorsf/reference/orsf_vint.md)
  to compute variable interaction scores using partial dependence.

- Added
  [`orsf_update()`](https://bcjaeger.github.io/aorsf/reference/orsf_update.md),
  which can copy and modify an `obliqueForest` or modify it in place.

- Added `orsf_control` functions for classification, regression, and
  survival (<https://github.com/ropensci/aorsf/pull/25>).

- optimization implemented for matrix multiplication during prediction
  (<https://github.com/ropensci/aorsf/pull/20>)

## aorsf 0.1.1

CRAN release: 2023-10-26

- Fixed an uninitialized value for `pd_type`

- Fixed various issues related to memory leaks

## aorsf 0.1.0

CRAN release: 2023-10-13

- Re-worked internal C++ routines following the design of `ranger`.

- Re-worked how progress is printed to console when `verbose_progress`
  is `TRUE`, following the design of `ranger`. Messages now indicate the
  action being taken, the % complete, and the approximate time until
  finishing the action.

- Improved variable importance, following the design of `ranger`.
  Importance is now computed tree-by-tree instead of by aggregate.
  Additionally, mortality is the type of prediction used for importance
  with survival trees, since mortality does not depend on
  `pred_horizon`.

- Allowed multi-threading to be performed in
  [`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md),
  `predict.orsf_fit()`, and functions in the
  [`orsf_vi()`](https://bcjaeger.github.io/aorsf/reference/orsf_vi.md)
  and `orsf_pd()` family.

- Allowed sampling without replacement and sampling a specific fraction
  of observations in
  [`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md)

- Included Harrell’s C-statistic as an option for assessing goodness of
  splits while growing trees.

- Fixed an issue where an uninformative error message would occur when
  `pred_horizon` was \> max(time) for `orsf_summarize_uni`. Thanks to
  [@JyHao1](https://github.com/JyHao1) and
  [@DustinMLong](https://github.com/DustinMLong) for finding this!

## aorsf 0.0.7

CRAN release: 2023-01-12

- Additional changes in internal testing to avoid problems with ATLAS

## aorsf 0.0.6

CRAN release: 2023-01-06

- Minor fix for internal tests that were failing when run on ATLAS

## aorsf 0.0.5

CRAN release: 2022-12-14

- [`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) no
  longer throws errors or warnings when you try to give it a single
  predictor. A note was added to the documentation in the details of
  [`?orsf`](https://bcjaeger.github.io/aorsf/reference/orsf.md) that
  explains why using a single predictor with
  [`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) is
  somewhat useless. This was done to resolve
  <https://github.com/mlr-org/mlr3extralearners/issues/259>.

- `predict.orsf_fit` now accepts `pred_horizon = 0` and returns sensible
  values. Thanks to [@mattwarkentin](https://github.com/mattwarkentin)
  for the feature request.

- added a function to perform variable selection,
  [`orsf_vs()`](https://bcjaeger.github.io/aorsf/reference/orsf_vs.md).

- Made variable importance consistent with respect to `group_factors`.
  Originally, the output from `orsf` would have ungrouped VI values
  while `orsf_vi` would have grouped values. With this update, `orsf`
  defaults to grouped values. The ungrouped values can still be
  recovered.

- Fixed an issue in `orsf_pd` functions where output data were not being
  returned on the original scale.

## aorsf 0.0.4

CRAN release: 2022-11-07

- `orsf` formulas now accepts `Surv` objects (see
  <https://github.com/ropensci/aorsf/issues/11>)

- Added `verbose_progress` input to `orsf`, which prints messages to
  console indicating progress.

- Allowance of missing values for `orsf`. Mean and mode imputation is
  performed for observations with missing data. These values can also be
  used to impute new data with missing values.

- Centering and scaling of predictors is now done prior to growing the
  forest.

## aorsf 0.0.3

CRAN release: 2022-10-09

- Included rOpenSci reviewers Christopher Jackson, Marvin N Wright, and
  Lukas Burk in `DESCRIPTION` as reviewers. Thank you!

- Added clarification to docs about pros/cons of different variable
  importance techniques

- Added regression tests for `aorsf` versus `obliqueRSF` (they should be
  similar)

- Additional support and tests for functions with long right hand sides

- Updated out-of-bag vignette with more appropriate custom functions.

- Allow status values in input data to be more general, i.e., not just 0
  and 1.

- Allow missing values in `predict` functions, including partial
  dependence.

## aorsf 0.0.2

CRAN release: 2022-09-05

- Modified unit tests for compatibility with extra checks run through
  CRAN.

## aorsf 0.0.1

CRAN release: 2022-08-23

- Added
  [`orsf_control_custom()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_custom.md),
  which allows users to submit custom functions for identifying linear
  combinations of inputs while growing oblique decision trees.

- Added `weights` input to `orsf`, allowing users to over or under fit
  `orsf` to specific data in their training set.

- Added `chf` and `mort` options to `predict.orsf_fit()`. Mortality
  predictions are not fully implemented yet - they are not supported in
  partial dependence or out-of-bag error estimates. These features will
  be added in a future update.

## aorsf 0.0.0.9000

- Core features implemented: fit, interpret, and predict using oblique
  random survival forests.

- Vignettes + Readme covering usage of core features.

- Website hosted through GitHub pages, managed with `pkgdown`.
