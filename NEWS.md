# aorsf 0.0.0.9000

* Core features implemented: fit, interpret, and predict using oblique random survival forests.

* Vignettes + Readme covering usage of core features.

* Website hosted through GitHub pages, managed with `pkgdown`.

# aorsf 0.0.1

* Added `orsf_control_custom()`, which allows users to submit custom functions for identifying linear combinations of inputs while growing oblique decision trees.

* Added `weights` input to `orsf`, allowing users to over or under fit `orsf` to specific data in their training set.

* Added `chf` and `mort` options to `predict.orsf_fit()`. Mortality predictions are not fully implemented yet - they are not supported in partial dependence or out-of-bag error estimates. These features will be added in a future update.

# aorsf 0.0.2

* Modified unit tests for compatibility with extra checks run through CRAN.
