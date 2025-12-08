# Prediction for ObliqueForest Objects

Compute predicted values from an oblique random forest. Predictions may
be returned in aggregate (i.e., averaging over all the trees) or
tree-specific.

## Usage

``` r
# S3 method for class 'ObliqueForest'
predict(
  object,
  new_data = NULL,
  pred_type = NULL,
  pred_horizon = NULL,
  pred_aggregate = TRUE,
  pred_simplify = FALSE,
  oobag = FALSE,
  na_action = NULL,
  boundary_checks = TRUE,
  n_thread = NULL,
  verbose_progress = NULL,
  ...
)
```

## Arguments

- object:

  (*ObliqueForest*) a trained oblique random forest object (see
  [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md)).

- new_data:

  a [data.frame](https://rdrr.io/r/base/data.frame.html),
  [tibble](https://tibble.tidyverse.org/reference/tibble-package.html),
  or
  [data.table](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
  to compute predictions in.

- pred_type:

  (*character*) the type of predictions to compute. Valid options for
  survival are:

  - 'risk' : probability of having an event at or before `pred_horizon`.

  - 'surv' : 1 - risk.

  - 'chf': cumulative hazard function

  - 'mort': mortality prediction

  - 'time': survival time prediction

  For classification:

  - 'prob': probability for each class

  - 'class': predicted class

  For regression:

  - 'mean': predicted mean, i.e., the expected value

- pred_horizon:

  (*double*) Only relevent for survival forests. A value or vector
  indicating the time(s) that predictions will be calibrated to. E.g.,
  if you were predicting risk of incident heart failure within the next
  10 years, then `pred_horizon = 10`. `pred_horizon` can be `NULL` if
  `pred_type` is `'mort'`, since mortality predictions are aggregated
  over all event times

- pred_aggregate:

  (*logical*) If `TRUE` (the default), predictions will be aggregated
  over all trees by taking the mean. If `FALSE`, the returned output
  will contain one row per observation and one column for each tree. If
  the length of `pred_horizon` is two or more and `pred_aggregate` is
  `FALSE`, then the result will be a list of such matrices, with the
  i'th item in the list corresponding to the i'th value of
  `pred_horizon`.

- pred_simplify:

  (*logical*) If `FALSE` (the default), predictions will always be
  returned in a numeric matrix or a list of numeric matrices. If `TRUE`,
  predictions may be simplified to a vector, e.g., if `pred_type` is
  `'mort'` for survival or `'class'` for classification, or an array of
  matrices if `length(pred_horizon) > 1`.

- oobag:

  (*logical*) If `FALSE` (the default), predictions will be computed
  using all trees for each observation. If `TRUE`, then out-of-bag
  predictions will be computed. This input parameter should only be set
  to `TRUE` if `new_data` is `NULL`.

- na_action:

  (*character*) what should happen when `new_data` contains missing
  values (i.e., `NA` values). Valid options are:

  - 'fail' : an error is thrown if `new_data` contains `NA` values

  - 'pass' : the output will have `NA` in all rows where `new_data` has
    1 or more `NA` value for the predictors used by `object`

  - 'omit' : rows in `new_data` with incomplete data will be dropped

  - 'impute_meanmode' : missing values for continuous and categorical
    variables in `new_data` will be imputed using the mean and mode,
    respectively. To clarify, the mean and mode used to impute missing
    values are from the training data of `object`, not from `new_data`.

- boundary_checks:

  (*logical*) if `TRUE`, `pred_horizon` will be checked to make sure the
  requested values are less than the maximum observed time in `object`'s
  training data. If `FALSE`, these checks are skipped.

- n_thread:

  (*integer*) number of threads to use while computing predictions.
  Default is 0, which allows a suitable number of threads to be used
  based on availability.

- verbose_progress:

  (*logical*) if `TRUE`, progress messages are printed in the console.
  If `FALSE` (the default), nothing is printed.

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

a `matrix` of predictions. Column `j` of the matrix corresponds to value
`j` in `pred_horizon`. Row `i` of the matrix corresponds to row `i` in
`new_data`.

## Details

`new_data` must have the same columns with equivalent types as the data
used to train `object`. Also, factors in `new_data` must not have levels
that were not in the data used to train `object`.

`pred_horizon` values should not exceed the maximum follow-up time in
`object`'s training data, but if you truly want to do this, set
`boundary_checks = FALSE` and you can use a `pred_horizon` as large as
you want. Note that predictions beyond the maximum follow-up time in the
`object`'s training data are equal to predictions at the maximum
follow-up time, because `aorsf` does not estimate survival beyond its
maximum observed time.

If unspecified, `pred_horizon` may be automatically specified as the
value used for `oobag_pred_horizon` when `object` was created (see
[orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md)).

## Examples

    library(aorsf)

### Classification

    set.seed(329)

    index_train <- sample(nrow(penguins_orsf), 150)

    penguins_orsf_train <- penguins_orsf[index_train, ]
    penguins_orsf_test <- penguins_orsf[-index_train, ]

    fit_clsf <- orsf(data = penguins_orsf_train,
                     formula = species ~ .)

Predict probability for each class or the predicted class:

    # predicted probabilities, the default
    predict(fit_clsf,
            new_data = penguins_orsf_test[1:5, ],
            pred_type = 'prob')

    ##         Adelie  Chinstrap      Gentoo
    ## [1,] 0.9405286 0.04125900 0.018212368
    ## [2,] 0.9628964 0.03459853 0.002505059
    ## [3,] 0.9029383 0.08527806 0.011783605
    ## [4,] 0.9301983 0.05180907 0.017992625
    ## [5,] 0.7968234 0.16538539 0.037791201

    # predicted class (as a matrix by default)
    predict(fit_clsf,
            new_data = penguins_orsf_test[1:5, ],
            pred_type = 'class')

    ##      [,1]
    ## [1,]    1
    ## [2,]    1
    ## [3,]    1
    ## [4,]    1
    ## [5,]    1

    # predicted class (as a factor if you use simplify)
    predict(fit_clsf,
            new_data = penguins_orsf_test[1:5, ],
            pred_type = 'class',
            pred_simplify = TRUE)

    ## [1] Adelie Adelie Adelie Adelie Adelie
    ## Levels: Adelie Chinstrap Gentoo

### Regression

    set.seed(329)

    index_train <- sample(nrow(penguins_orsf), 150)

    penguins_orsf_train <- penguins_orsf[index_train, ]
    penguins_orsf_test <- penguins_orsf[-index_train, ]

    fit_regr <- orsf(data = penguins_orsf_train,
                     formula = bill_length_mm ~ .)

Predict the mean value of the outcome:

    predict(fit_regr,
            new_data = penguins_orsf_test[1:5, ],
            pred_type = 'mean')

    ##          [,1]
    ## [1,] 37.74136
    ## [2,] 37.42367
    ## [3,] 37.04598
    ## [4,] 39.89602
    ## [5,] 39.14848

### Survival

Begin by fitting an oblique survival random forest:

    set.seed(329)

    index_train <- sample(nrow(pbc_orsf), 150)

    pbc_orsf_train <- pbc_orsf[index_train, ]
    pbc_orsf_test <- pbc_orsf[-index_train, ]

    fit_surv <- orsf(data = pbc_orsf_train,
                     formula = Surv(time, status) ~ . - id,
                     oobag_pred_horizon = 365.25 * 5)

Predict risk, survival, or cumulative hazard at one or several times:

    # predicted risk, the default
    predict(fit_surv,
            new_data = pbc_orsf_test[1:5, ],
            pred_type = 'risk',
            pred_horizon = c(500, 1000, 1500))

    ##             [,1]        [,2]       [,3]
    ## [1,] 0.013648562 0.058393393 0.11184029
    ## [2,] 0.003811413 0.026857586 0.04774151
    ## [3,] 0.030548361 0.100600301 0.14847107
    ## [4,] 0.040381075 0.169596943 0.27018952
    ## [5,] 0.001484698 0.006663576 0.01337655

    # predicted survival, i.e., 1 - risk
    predict(fit_surv,
            new_data = pbc_orsf_test[1:5, ],
            pred_type = 'surv',
            pred_horizon = c(500, 1000, 1500))

    ##           [,1]      [,2]      [,3]
    ## [1,] 0.9863514 0.9416066 0.8881597
    ## [2,] 0.9961886 0.9731424 0.9522585
    ## [3,] 0.9694516 0.8993997 0.8515289
    ## [4,] 0.9596189 0.8304031 0.7298105
    ## [5,] 0.9985153 0.9933364 0.9866235

    # predicted cumulative hazard function
    # (expected number of events for person i at time j)
    predict(fit_surv,
            new_data = pbc_orsf_test[1:5, ],
            pred_type = 'chf',
            pred_horizon = c(500, 1000, 1500))

    ##             [,1]        [,2]       [,3]
    ## [1,] 0.015395388 0.067815817 0.14942956
    ## [2,] 0.004022524 0.028740305 0.05424314
    ## [3,] 0.034832754 0.127687156 0.20899732
    ## [4,] 0.059978334 0.233048809 0.42562310
    ## [5,] 0.001651365 0.007173177 0.01393016

Predict mortality, defined as the number of events in the forest’s
population if all observations had characteristics like the current
observation. This type of prediction does not require you to specify a
prediction horizon

    predict(fit_surv,
            new_data = pbc_orsf_test[1:5, ],
            pred_type = 'mort')

    ##           [,1]
    ## [1,] 23.405016
    ## [2,] 15.362916
    ## [3,] 26.180648
    ## [4,] 36.515629
    ## [5,]  5.856674
