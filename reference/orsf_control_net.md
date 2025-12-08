# Penalized Cox regression ORSF control

Use regularized Cox proportional hazard models to identify linear
combinations of input variables while fitting an
[orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md) model.

## Usage

``` r
orsf_control_net(alpha = 1/2, df_target = NULL, ...)
```

## Arguments

- alpha:

  (*double*) The elastic net mixing parameter. A value of 1 gives the
  lasso penalty, and a value of 0 gives the ridge penalty. If multiple
  values of alpha are given, then a penalized model is fit using each
  alpha value prior to splitting a node.

- df_target:

  (*integer*) Preferred number of variables used in a linear
  combination.

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

an object of class `'orsf_control'`, which should be used as an input
for the `control` argument of
[orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md).

## Details

**\[superseded\]**

`df_target` has to be less than `mtry`, which is a separate argument in
[orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md) that
indicates the number of variables chosen at random prior to finding a
linear combination of those variables.

## References

1.  Simon, Noah, Friedman, Jerome, Hastie, Trevor, Tibshirani, Rob
    (2011). "Regularization paths for Cox's proportional hazards model
    via coordinate descent." *Journal of statistical software*, *39*(5),
    1.

## See also

linear combination control functions
[`orsf_control()`](https://bcjaeger.github.io/aorsf/reference/orsf_control.md),
[`orsf_control_cph()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_cph.md),
[`orsf_control_custom()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_custom.md),
[`orsf_control_fast()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_fast.md)
