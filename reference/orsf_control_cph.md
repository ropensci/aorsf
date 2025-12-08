# Cox regression ORSF control

Use the coefficients from a proportional hazards model to create linear
combinations of predictor variables while fitting an
[orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md) model.

## Usage

``` r
orsf_control_cph(method = "efron", eps = 1e-09, iter_max = 20, ...)
```

## Arguments

- method:

  (*character*) a character string specifying the method for tie
  handling. If there are no ties, all the methods are equivalent. Valid
  options are 'breslow' and 'efron'. The Efron approximation is the
  default because it is more accurate when dealing with tied event times
  and has similar computational efficiency compared to the Breslow
  method.

- eps:

  (*double*) When using Newton Raphson scoring to identify linear
  combinations of inputs, iteration continues in the algorithm until the
  relative change in the log partial likelihood is less than `eps`, or
  the absolute change is less than `sqrt(eps)`. Must be positive. A
  default value of 1e-09 is used for consistency with
  [survival::coxph.control](https://rdrr.io/pkg/survival/man/coxph.control.html).

- iter_max:

  (*integer*) iteration continues until convergence (see `eps` above) or
  the number of attempted iterations is equal to `iter_max`.

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

an object of class `'orsf_control'`, which should be used as an input
for the `control` argument of
[orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md).

## Details

**\[superseded\]**

code from the [survival
package](https://github.com/therneau/survival/blob/master/src/coxfit6.c)
was modified to make this routine.

For more details on the Cox proportional hazards model, see
[coxph](https://rdrr.io/pkg/survival/man/coxph.html) and/or Therneau and
Grambsch (2000).

## References

Therneau T.M., Grambsch P.M. (2000) The Cox Model. In: Modeling Survival
Data: Extending the Cox Model. Statistics for Biology and Health.
Springer, New York, NY. DOI: 10.1007/978-1-4757-3294-8_3

## See also

linear combination control functions
[`orsf_control()`](https://bcjaeger.github.io/aorsf/reference/orsf_control.md),
[`orsf_control_custom()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_custom.md),
[`orsf_control_fast()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_fast.md),
[`orsf_control_net()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_net.md)
