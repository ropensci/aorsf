# Variable selection

Variable selection

## Usage

``` r
orsf_vs(
  object,
  n_predictor_min = 3,
  n_predictor_drop = 1,
  verbose_progress = NULL
)
```

## Arguments

- object:

  (*ObliqueForest*) a trained oblique random forest object (see
  [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md)).

- n_predictor_min:

  (*integer*) the minimum number of predictors allowed

- n_predictor_drop:

  (*integer*) the number of predictors dropped at each step

- verbose_progress:

  (*logical*) not implemented yet. Should progress be printed to the
  console?

## Value

a [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) with
four columns:

- *n_predictors*: the number of predictors used

- *stat_value*: the out-of-bag statistic

- *variables_included*: the names of the variables included

- *predictors_included*: the names of the predictors included

- *predictor_dropped*: the predictor selected to be dropped

## Details

The difference between `variables_included` and `predictors_included` is
referent coding. The `variable` would be the name of a factor variable
in the training data, while the `predictor` would be the name of that
same factor with the levels of the factor appended. For example, if the
variable is `diabetes` with `levels = c("no", "yes")`, then the variable
name is `diabetes` and the predictor name is `diabetes_yes`.

`tree_seeds` should be specified in `object` so that each successive run
of `orsf` will be evaluated in the same out-of-bag samples as the
initial run.

## Examples

``` r
object <- orsf(formula = time + status ~ .,
               data = pbc_orsf,
               n_tree = 25,
               importance = 'anova')

orsf_vs(object, n_predictor_min = 15)
#>    n_predictors stat_value                            variables_included
#>           <int>      <num>                                        <list>
#> 1:           15  0.8244179     age,albumin,ascites,ast,bili,chol,...[14]
#> 2:           16  0.8395229 age,albumin,alk.phos,ascites,ast,bili,...[15]
#> 3:           17  0.8429085 age,albumin,alk.phos,ascites,ast,bili,...[16]
#> 4:           18  0.8340539 age,albumin,alk.phos,ascites,ast,bili,...[17]
#> 5:           19  0.8210323 age,albumin,alk.phos,ascites,ast,bili,...[18]
#>                                      predictors_included predictor_dropped
#>                                                   <list>            <list>
#> 1: id,age,ascites_1,hepato_1,spiders_1,edema_0.5,...[15]                NA
#> 2: id,age,ascites_1,hepato_1,spiders_1,edema_0.5,...[16]          alk.phos
#> 3:     id,age,sex_f,ascites_1,hepato_1,spiders_1,...[17]             sex_f
#> 4:     id,age,sex_f,ascites_1,hepato_1,spiders_1,...[18]          platelet
#> 5:   id,trt_placebo,age,sex_f,ascites_1,hepato_1,...[19]       trt_placebo
```
