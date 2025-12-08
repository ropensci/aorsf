# Package index

## Oblique random forests (RFs)

Fit, inspect, summarize, and apply oblique RFs

- [`orsf()`](https://bcjaeger.github.io/aorsf/reference/orsf.md)
  [`orsf_train()`](https://bcjaeger.github.io/aorsf/reference/orsf.md) :
  Oblique Random Forests
- [`orsf_update()`](https://bcjaeger.github.io/aorsf/reference/orsf_update.md)
  : Update Forest Parameters
- [`print(`*`<ObliqueForest>`*`)`](https://bcjaeger.github.io/aorsf/reference/print.ObliqueForest.md)
  : Inspect Forest Parameters
- [`predict(`*`<ObliqueForest>`*`)`](https://bcjaeger.github.io/aorsf/reference/predict.ObliqueForest.md)
  : Prediction for ObliqueForest Objects
- [`orsf_summarize_uni()`](https://bcjaeger.github.io/aorsf/reference/orsf_summarize_uni.md)
  : Univariate summary
- [`print(`*`<orsf_summary_uni>`*`)`](https://bcjaeger.github.io/aorsf/reference/print.orsf_summary_uni.md)
  : Print ORSF summary

## Control how your oblique RF works

Choose how to identify linear combinations of predictors and set tuning
parameters for your approach

- [`orsf_control()`](https://bcjaeger.github.io/aorsf/reference/orsf_control.md)
  [`orsf_control_classification()`](https://bcjaeger.github.io/aorsf/reference/orsf_control.md)
  [`orsf_control_regression()`](https://bcjaeger.github.io/aorsf/reference/orsf_control.md)
  [`orsf_control_survival()`](https://bcjaeger.github.io/aorsf/reference/orsf_control.md)
  : Oblique random forest control
- [`orsf_control_cph()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_cph.md)
  : Cox regression ORSF control
- [`orsf_control_custom()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_custom.md)
  **\[superseded\]** : Custom ORSF control
- [`orsf_control_fast()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_fast.md)
  : Accelerated ORSF control
- [`orsf_control_net()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_net.md)
  : Penalized Cox regression ORSF control

## Variable importance/selection

Estimate the importance of individual variables and conduct variable
selection using ORSFs

- [`orsf_vi()`](https://bcjaeger.github.io/aorsf/reference/orsf_vi.md)
  [`orsf_vi_negate()`](https://bcjaeger.github.io/aorsf/reference/orsf_vi.md)
  [`orsf_vi_permute()`](https://bcjaeger.github.io/aorsf/reference/orsf_vi.md)
  [`orsf_vi_anova()`](https://bcjaeger.github.io/aorsf/reference/orsf_vi.md)
  : Variable Importance
- [`orsf_vint()`](https://bcjaeger.github.io/aorsf/reference/orsf_vint.md)
  : Variable Interactions
- [`orsf_vs()`](https://bcjaeger.github.io/aorsf/reference/orsf_vs.md) :
  Variable selection

## Partial dependence and individual conditional expectations

Interpret your model by generating partial dependence or individual
conditional expectation values. Plotting functions not included (but see
examples)

- [`orsf_ice_oob()`](https://bcjaeger.github.io/aorsf/reference/orsf_ice_oob.md)
  [`orsf_ice_inb()`](https://bcjaeger.github.io/aorsf/reference/orsf_ice_oob.md)
  [`orsf_ice_new()`](https://bcjaeger.github.io/aorsf/reference/orsf_ice_oob.md)
  : Individual Conditional Expectations
- [`orsf_pd_oob()`](https://bcjaeger.github.io/aorsf/reference/orsf_pd_oob.md)
  [`orsf_pd_inb()`](https://bcjaeger.github.io/aorsf/reference/orsf_pd_oob.md)
  [`orsf_pd_new()`](https://bcjaeger.github.io/aorsf/reference/orsf_pd_oob.md)
  : Partial dependence
- [`pred_spec_auto()`](https://bcjaeger.github.io/aorsf/reference/pred_spec_auto.md)
  : Automatic variable values for dependence

## Example survival data

Datasets used in examples and vignettes.

- [`pbc_orsf`](https://bcjaeger.github.io/aorsf/reference/pbc_orsf.md) :
  Mayo Clinic Primary Biliary Cholangitis Data
- [`penguins_orsf`](https://bcjaeger.github.io/aorsf/reference/penguins_orsf.md)
  : Size measurements for adult foraging penguins near Palmer Station,
  Antarctica

## Miscellaneous

Functions that don’t fit neatly into a category above, but are still
helpful.

- [`as.data.table(`*`<orsf_summary_uni>`*`)`](https://bcjaeger.github.io/aorsf/reference/as.data.table.orsf_summary_uni.md)
  : Coerce to data.table
- [`orsf_time_to_train()`](https://bcjaeger.github.io/aorsf/reference/orsf_time_to_train.md)
  : Estimate training time

## Back-end functions

Techniques used by aorsf that may be helpful in other contexts.

- [`orsf_scale_cph()`](https://bcjaeger.github.io/aorsf/reference/orsf_scale_cph.md)
  [`orsf_unscale_cph()`](https://bcjaeger.github.io/aorsf/reference/orsf_scale_cph.md)
  : Scale input data
