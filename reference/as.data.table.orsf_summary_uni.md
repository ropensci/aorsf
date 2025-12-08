# Coerce to data.table

Convert an 'orsf_summary' object into a `data.table` object.

## Usage

``` r
# S3 method for class 'orsf_summary_uni'
as.data.table(x, ...)
```

## Arguments

- x:

  an object of class 'orsf_summary_uni'

- ...:

  not used

## Value

a
[data.table](https://rdatatable.gitlab.io/data.table/reference/data.table.html)

## Examples

``` r
if (FALSE) { # \dontrun{

library(data.table)

object <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 25)

smry <- orsf_summarize_uni(object, n_variables = 2)

as.data.table(smry)

} # }

```
