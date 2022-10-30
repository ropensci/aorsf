# aorsf 0.0.4

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

