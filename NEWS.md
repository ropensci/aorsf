# aorsf 0.1.1.9000

* Added `orsf_control` functions for classification, regression, and survival.

* optimization implemented for matrix multiplication during prediction (https://github.com/ropensci/aorsf/pull/20)

# aorsf 0.1.1

* fixed an uninitialized value for `pd_type`

# aorsf 0.1.0

* Re-worked internal C++ routines following the design of `ranger`. 

* Re-worked how progress is printed to console when `verbose_progress` is `TRUE`, following the design of `ranger`. Messages now indicate the action being taken, the % complete, and the approximate time until finishing the action. 

* Improved variable importance, following the design of `ranger`. Importance is now computed tree-by-tree instead of by aggregate. Additionally, mortality is the type of prediction used for importance with survival trees, since mortality does not depend on `pred_horizon`.

* Allowed multi-threading to be performed in `orsf()`, `predict.orsf_fit()`, and functions in the `orsf_vi()` and `orsf_pd()` family.

* Allowed sampling without replacement and sampling a specific fraction of observations in `orsf()`

* Included Harrell's C-statistic as an option for assessing goodness of splits while growing trees.

* Fixed an issue where an uninformative error message would occur when `pred_horizon` was > max(time) for `orsf_summarize_uni`. Thanks to @JyHao1 and @DustinMLong for finding this!

# aorsf 0.0.7

* Additional changes in internal testing to avoid problems with ATLAS

# aorsf 0.0.6

* Minor fix for internal tests that were failing when run on ATLAS

# aorsf 0.0.5

* `orsf()` no longer throws errors or warnings when you try to give it a single predictor. A note was added to the documentation in the details of `?orsf` that explains why using a single predictor with `orsf()` is somewhat useless. This was done to resolve https://github.com/mlr-org/mlr3extralearners/issues/259.

* `predict.orsf_fit` now accepts `pred_horizon = 0` and returns sensible values. Thanks to @mattwarkentin for the feature request.

* added a function to perform variable selection, `orsf_vs()`.

* Made variable importance consistent with respect to `group_factors`. Originally, the output from `orsf` would have ungrouped VI values while `orsf_vi` would have grouped values. With this update, `orsf` defaults to grouped values. The ungrouped values can still be recovered.

* Fixed an issue in `orsf_pd` functions where output data were not being returned on the original scale.


# aorsf 0.0.4

* `orsf` formulas now accepts `Surv` objects (see https://github.com/ropensci/aorsf/issues/11)

* Added `verbose_progress` input to `orsf`, which prints messages to console indicating progress. 

* Allowance of missing values for `orsf`. Mean and mode imputation is performed for observations with missing data. These values can also be used to impute new data with missing values.

* Centering and scaling of predictors is now done prior to growing the forest. 

# aorsf 0.0.3

* Included rOpenSci reviewers Christopher Jackson, Marvin N Wright, and Lukas Burk in `DESCRIPTION` as reviewers. Thank you!

* Added clarification to docs about pros/cons of different variable importance techniques

* Added regression tests for `aorsf` versus `obliqueRSF` (they should be similar)

* Additional support and tests for functions with long right hand sides

* Updated out-of-bag vignette with more appropriate custom functions.

* Allow status values in input data to be more general, i.e., not just 0 and 1.

* Allow missing values in `predict` functions, including partial dependence.

# aorsf 0.0.2

* Modified unit tests for compatibility with extra checks run through CRAN.

# aorsf 0.0.1

* Added `orsf_control_custom()`, which allows users to submit custom functions for identifying linear combinations of inputs while growing oblique decision trees.

* Added `weights` input to `orsf`, allowing users to over or under fit `orsf` to specific data in their training set.

* Added `chf` and `mort` options to `predict.orsf_fit()`. Mortality predictions are not fully implemented yet - they are not supported in partial dependence or out-of-bag error estimates. These features will be added in a future update.

# aorsf 0.0.0.9000

* Core features implemented: fit, interpret, and predict using oblique random survival forests.

* Vignettes + Readme covering usage of core features.

* Website hosted through GitHub pages, managed with `pkgdown`.

