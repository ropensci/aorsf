# Custom ORSF control

**\[superseded\]**

## Usage

``` r
orsf_control_custom(beta_fun, ...)
```

## Arguments

- beta_fun:

  (*function*) a function to define coefficients used in linear
  combinations of predictor variables. `beta_fun` must accept three
  inputs named `x_node`, `y_node` and `w_node`, and should expect the
  following types and dimensions:

  - `x_node` (*matrix*; `n` rows, `p` columns)

  - `y_node` (*matrix*; `n` rows, `2` columns)

  - `w_node` (*matrix*; `n` rows, `1` column)

  In addition, `beta_fun` must return a matrix with p rows and 1 column.
  If any of these conditions are not met, `orsf_control_custom()` will
  let you know.

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

an object of class `'orsf_control'`, which should be used as an input
for the `control` argument of
[orsf](https://bcjaeger.github.io/aorsf/reference/orsf.md).

## See also

linear combination control functions
[`orsf_control()`](https://bcjaeger.github.io/aorsf/reference/orsf_control.md),
[`orsf_control_cph()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_cph.md),
[`orsf_control_fast()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_fast.md),
[`orsf_control_net()`](https://bcjaeger.github.io/aorsf/reference/orsf_control_net.md)
