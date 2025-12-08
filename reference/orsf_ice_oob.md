# Individual Conditional Expectations

Compute individual conditional expectations for an oblique random
forest. Unlike partial dependence, which shows the expected prediction
as a function of one or multiple predictors, individual conditional
expectations (ICE) show the prediction for an individual observation as
a function of a predictor. You can compute individual conditional
expectations three ways using a random forest:

- using in-bag predictions for the training data

- using out-of-bag predictions for the training data

- using predictions for a new set of data

See examples for more details

## Usage

``` r
orsf_ice_oob(
  object,
  pred_spec,
  pred_horizon = NULL,
  pred_type = NULL,
  expand_grid = TRUE,
  boundary_checks = TRUE,
  n_thread = NULL,
  verbose_progress = NULL,
  ...
)

orsf_ice_inb(
  object,
  pred_spec,
  pred_horizon = NULL,
  pred_type = NULL,
  expand_grid = TRUE,
  boundary_checks = TRUE,
  n_thread = NULL,
  verbose_progress = NULL,
  ...
)

orsf_ice_new(
  object,
  pred_spec,
  new_data,
  pred_horizon = NULL,
  pred_type = NULL,
  na_action = "fail",
  expand_grid = TRUE,
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

- pred_spec:

  (*named list*, *pspec_auto*, or *data.frame*).

  - If `pred_spec` is a named list, Each item in the list should be a
    vector of values that will be used as points in the partial
    dependence function. The name of each item in the list should
    indicate which variable will be modified to take the corresponding
    values.

  - If `pred_spec` is created using
    [`pred_spec_auto()`](https://bcjaeger.github.io/aorsf/reference/pred_spec_auto.md),
    all that is needed is the names of variables to use (see
    [pred_spec_auto](https://bcjaeger.github.io/aorsf/reference/pred_spec_auto.md)).

  - If `pred_spec` is a `data.frame`, columns will indicate variable
    names, values will indicate variable values, and partial dependence
    will be computed using the inputs on each row.

- pred_horizon:

  (*double*) Only relevent for survival forests. A value or vector
  indicating the time(s) that predictions will be calibrated to. E.g.,
  if you were predicting risk of incident heart failure within the next
  10 years, then `pred_horizon = 10`. `pred_horizon` can be `NULL` if
  `pred_type` is `'mort'`, since mortality predictions are aggregated
  over all event times

- pred_type:

  (*character*) the type of predictions to compute. Valid Valid options
  for survival are:

  - 'risk' : probability of having an event at or before `pred_horizon`.

  - 'surv' : 1 - risk.

  - 'chf': cumulative hazard function

  - 'mort': mortality prediction

  - 'time': survival time prediction

  For classification:

  - 'prob': probability for each class

  For regression:

  - 'mean': predicted mean, i.e., the expected value

- expand_grid:

  (*logical*) if `TRUE`, partial dependence will be computed at all
  possible combinations of inputs in `pred_spec`. If `FALSE`, partial
  dependence will be computed for each variable in `pred_spec`,
  separately.

- boundary_checks:

  (*logical*) if `TRUE`, `pred_spec` will be checked to make sure the
  requested values are between the 10th and 90th percentile in the
  object's training data. If `FALSE`, these checks are skipped.

- n_thread:

  (*integer*) number of threads to use while computing predictions.
  Default is 0, which allows a suitable number of threads to be used
  based on availability.

- verbose_progress:

  (*logical*) if `TRUE`, progress will be printed to console. If `FALSE`
  (the default), nothing will be printed.

- ...:

  Further arguments passed to or from other methods (not currently
  used).

- new_data:

  a [data.frame](https://rdrr.io/r/base/data.frame.html),
  [tibble](https://tibble.tidyverse.org/reference/tibble-package.html),
  or
  [data.table](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
  to compute predictions in.

- na_action:

  (*character*) what should happen when `new_data` contains missing
  values (i.e., `NA` values). Valid options are:

  - 'fail' : an error is thrown if `new_data` contains `NA` values

  - 'omit' : rows in `new_data` with incomplete data will be dropped

## Value

a
[data.table](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
containing individual conditional expectations for the specified
variable(s) and, if relevant, at the specified prediction horizon(s).

## Examples

You can compute individual conditional expectation and individual
conditional expectations in three ways:

- using in-bag predictions for the training data. In-bag individual
  conditional expectation indicates relationships that the model has
  learned during training. This is helpful if your goal is to interpret
  the model.

- using out-of-bag predictions for the training data. Out-of-bag
  individual conditional expectation indicates relationships that the
  model has learned during training but using the out-of-bag data
  simulates application of the model to new data. This is helpful if you
  want to test your model’s reliability or fairness in new data but you
  don’t have access to a large testing set.

- using predictions for a new set of data. New data individual
  conditional expectation shows how the model predicts outcomes for
  observations it has not seen. This is helpful if you want to test your
  model’s reliability or fairness.

### Classification

Begin by fitting an oblique classification random forest:

    set.seed(329)

    index_train <- sample(nrow(penguins_orsf), 150)

    penguins_orsf_train <- penguins_orsf[index_train, ]
    penguins_orsf_test <- penguins_orsf[-index_train, ]

    fit_clsf <- orsf(data = penguins_orsf_train,
                     formula = species ~ .)

Compute individual conditional expectation using out-of-bag data for
`flipper_length_mm = c(190, 210)`.

    pred_spec <- list(flipper_length_mm = c(190, 210))

    ice_oob <- orsf_ice_oob(fit_clsf, pred_spec = pred_spec)

    ice_oob

    ## Key: <class>
    ##      id_variable id_row  class flipper_length_mm       pred
    ##            <int> <char> <fctr>             <num>      <num>
    ##   1:           1      1 Adelie               190 0.92059968
    ##   2:           1      2 Adelie               190 0.80953569
    ##   3:           1      3 Adelie               190 0.84869374
    ##   4:           1      4 Adelie               190 0.93559660
    ##   5:           1      5 Adelie               190 0.97708693
    ##  ---
    ## 896:           2    146 Gentoo               210 0.25636964
    ## 897:           2    147 Gentoo               210 0.04798334
    ## 898:           2    148 Gentoo               210 0.07945140
    ## 899:           2    149 Gentoo               210 0.84811899
    ## 900:           2    150 Gentoo               210 0.10695367

There are two identifiers in the output:

- `id_variable` is an identifier for the current value of the
  variable(s) that are in the data. It is redundant if you only have one
  variable, but helpful if there are multiple variables.

- `id_row` is an identifier for the observation in the original data.

Note that predicted probabilities are returned for each class and each
observation in the data. Predicted probabilities for a given observation
and given variable value sum to 1. For example,

    ice_oob %>%
     .[flipper_length_mm == 190] %>%
     .[id_row == 1] %>%
     .[['pred']] %>%
     sum()

    ## [1] 1

### Regression

Begin by fitting an oblique regression random forest:

    set.seed(329)

    index_train <- sample(nrow(penguins_orsf), 150)

    penguins_orsf_train <- penguins_orsf[index_train, ]
    penguins_orsf_test <- penguins_orsf[-index_train, ]

    fit_regr <- orsf(data = penguins_orsf_train,
                     formula = bill_length_mm ~ .)

Compute individual conditional expectation using new data for
`flipper_length_mm = c(190, 210)`.

    pred_spec <- list(flipper_length_mm = c(190, 210))

    ice_new <- orsf_ice_new(fit_regr,
                            pred_spec = pred_spec,
                            new_data = penguins_orsf_test)

    ice_new

    ##      id_variable id_row flipper_length_mm     pred
    ##            <int> <char>             <num>    <num>
    ##   1:           1      1               190 37.94483
    ##   2:           1      2               190 37.61595
    ##   3:           1      3               190 37.53681
    ##   4:           1      4               190 39.49476
    ##   5:           1      5               190 38.95635
    ##  ---
    ## 362:           2    179               210 51.80471
    ## 363:           2    180               210 47.27183
    ## 364:           2    181               210 47.05031
    ## 365:           2    182               210 50.39028
    ## 366:           2    183               210 48.44774

You can also let `pred_spec_auto` pick reasonable values like so:

    pred_spec = pred_spec_auto(species, island, body_mass_g)

    ice_new <- orsf_ice_new(fit_regr,
                            pred_spec = pred_spec,
                            new_data = penguins_orsf_test)

    ice_new

    ##       id_variable id_row species    island body_mass_g     pred
    ##             <int> <char>  <fctr>    <fctr>       <num>    <num>
    ##    1:           1      1  Adelie    Biscoe        3200 37.78339
    ##    2:           1      2  Adelie    Biscoe        3200 37.73273
    ##    3:           1      3  Adelie    Biscoe        3200 37.71248
    ##    4:           1      4  Adelie    Biscoe        3200 40.25782
    ##    5:           1      5  Adelie    Biscoe        3200 40.04074
    ##   ---
    ## 8231:          45    179  Gentoo Torgersen        5300 46.14559
    ## 8232:          45    180  Gentoo Torgersen        5300 43.98050
    ## 8233:          45    181  Gentoo Torgersen        5300 44.59837
    ## 8234:          45    182  Gentoo Torgersen        5300 44.85146
    ## 8235:          45    183  Gentoo Torgersen        5300 44.23710

By default, all combinations of all variables are used. However, you can
also look at the variables one by one, separately, like so:

    ice_new <- orsf_ice_new(fit_regr,
                            expand_grid = FALSE,
                            pred_spec = pred_spec,
                            new_data = penguins_orsf_test)

    ice_new

    ##       id_variable id_row    variable value  level     pred
    ##             <int> <char>      <char> <num> <char>    <num>
    ##    1:           1      1     species    NA Adelie 37.74136
    ##    2:           1      2     species    NA Adelie 37.42367
    ##    3:           1      3     species    NA Adelie 37.04598
    ##    4:           1      4     species    NA Adelie 39.89602
    ##    5:           1      5     species    NA Adelie 39.14848
    ##   ---
    ## 2009:           5    179 body_mass_g  5300   <NA> 51.50196
    ## 2010:           5    180 body_mass_g  5300   <NA> 47.27055
    ## 2011:           5    181 body_mass_g  5300   <NA> 48.34064
    ## 2012:           5    182 body_mass_g  5300   <NA> 48.75828
    ## 2013:           5    183 body_mass_g  5300   <NA> 48.11020

And you can also bypass all the bells and whistles by using your own
`data.frame` for a `pred_spec`. (Just make sure you request values that
exist in the training data.)

    custom_pred_spec <- data.frame(species = 'Adelie',
                                   island = 'Biscoe')

    ice_new <- orsf_ice_new(fit_regr,
                            pred_spec = custom_pred_spec,
                            new_data = penguins_orsf_test)

    ice_new

    ##      id_variable id_row species island     pred
    ##            <int> <char>  <fctr> <fctr>    <num>
    ##   1:           1      1  Adelie Biscoe 38.52327
    ##   2:           1      2  Adelie Biscoe 38.32073
    ##   3:           1      3  Adelie Biscoe 37.71248
    ##   4:           1      4  Adelie Biscoe 41.68380
    ##   5:           1      5  Adelie Biscoe 40.91140
    ##  ---
    ## 179:           1    179  Adelie Biscoe 43.09493
    ## 180:           1    180  Adelie Biscoe 38.79455
    ## 181:           1    181  Adelie Biscoe 39.37734
    ## 182:           1    182  Adelie Biscoe 40.71952
    ## 183:           1    183  Adelie Biscoe 39.34501

### Survival

Begin by fitting an oblique survival random forest:

    set.seed(329)

    index_train <- sample(nrow(pbc_orsf), 150)

    pbc_orsf_train <- pbc_orsf[index_train, ]
    pbc_orsf_test <- pbc_orsf[-index_train, ]

    fit_surv <- orsf(data = pbc_orsf_train,
                     formula = Surv(time, status) ~ . - id,
                     oobag_pred_horizon = 365.25 * 5)

Compute individual conditional expectation using in-bag data for
`bili = c(1,2,3,4,5)`:

    ice_train <- orsf_ice_inb(fit_surv, pred_spec = list(bili = 1:5))
    ice_train

    ##      id_variable id_row pred_horizon  bili      pred
    ##            <int> <char>        <num> <num>     <num>
    ##   1:           1      1      1826.25     1 0.1290317
    ##   2:           1      2      1826.25     1 0.1242352
    ##   3:           1      3      1826.25     1 0.0963452
    ##   4:           1      4      1826.25     1 0.1172367
    ##   5:           1      5      1826.25     1 0.2030256
    ##  ---
    ## 746:           5    146      1826.25     5 0.7868537
    ## 747:           5    147      1826.25     5 0.2012954
    ## 748:           5    148      1826.25     5 0.4893605
    ## 749:           5    149      1826.25     5 0.4698220
    ## 750:           5    150      1826.25     5 0.9557285

If you don’t have specific values of a variable in mind, let
`pred_spec_auto` pick for you:

    ice_train <- orsf_ice_inb(fit_surv, pred_spec_auto(bili))
    ice_train

    ##      id_variable id_row pred_horizon  bili       pred
    ##            <int> <char>        <num> <num>      <num>
    ##   1:           1      1      1826.25  0.59 0.11706741
    ##   2:           1      2      1826.25  0.59 0.11562173
    ##   3:           1      3      1826.25  0.59 0.09110739
    ##   4:           1      4      1826.25  0.59 0.10069721
    ##   5:           1      5      1826.25  0.59 0.18769751
    ##  ---
    ## 746:           5    146      1826.25  7.21 0.82600898
    ## 747:           5    147      1826.25  7.21 0.29156437
    ## 748:           5    148      1826.25  7.21 0.58220919
    ## 749:           5    149      1826.25  7.21 0.54168688
    ## 750:           5    150      1826.25  7.21 0.96204106

Specify `pred_horizon` to get individual conditional expectation at each
value:

    ice_train <- orsf_ice_inb(fit_surv, pred_spec_auto(bili),
                              pred_horizon = seq(500, 3000, by = 500))
    ice_train

    ##       id_variable id_row pred_horizon  bili        pred
    ##             <int> <char>        <num> <num>       <num>
    ##    1:           1      1          500  0.59 0.008276627
    ##    2:           1      1         1000  0.59 0.055715858
    ##    3:           1      1         1500  0.59 0.084987224
    ##    4:           1      1         2000  0.59 0.123090885
    ##    5:           1      1         2500  0.59 0.165214938
    ##   ---
    ## 4496:           5    150         1000  7.21 0.835895969
    ## 4497:           5    150         1500  7.21 0.932657591
    ## 4498:           5    150         2000  7.21 0.965944498
    ## 4499:           5    150         2500  7.21 0.970325309
    ## 4500:           5    150         3000  7.21 0.979051377

Multi-prediction horizon ice comes with minimal extra computational
cost. Use a fine grid of time values and assess whether predictors have
time-varying effects.
