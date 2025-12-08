# Accelerated ORSF control

Fast methods to identify linear combinations of predictors while fitting
an [orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md) model.

## Usage

``` r
orsf_control_fast(method = "efron", do_scale = TRUE, ...)
```

## Arguments

- method:

  (*character*) a character string specifying the method for tie
  handling. If there are no ties, all the methods are equivalent. Valid
  options are 'breslow' and 'efron'. The Efron approximation is the
  default because it is more accurate when dealing with tied event times
  and has similar computational efficiency compared to the Breslow
  method.

- do_scale:

  (*logical*) if `TRUE`, values of predictors will be scaled prior to
  each instance of Newton Raphson scoring, using summary values from the
  data in the current node of the decision tree.

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

an object of class `'orsf_control'`, which should be used as an input
for the `control` argument of
[orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md).

## Details

code from the [survival
package](https://github.com/therneau/survival/blob/master/src/coxfit6.c)
was modified to make this routine.

Adjust `do_scale` *at your own risk*. Setting `do_scale = FALSE` will
reduce computation time but will also make the `orsf` model dependent on
the scale of your data, which is why the default value is `TRUE`.

## See also

linear combination control functions
[`orsf_control()`](https://bcjaeger.github.io/aorsf/reference/orsf_control.md),
[`orsf_control_cph()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_cph.md),
[`orsf_control_custom()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_custom.md),
[`orsf_control_net()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_net.md)
