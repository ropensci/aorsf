# Mayo Clinic Primary Biliary Cholangitis Data

These data are a light modification of the
[survival::pbc](https://rdrr.io/pkg/survival/man/pbc.html) data. The
modifications are:

## Usage

``` r
pbc_orsf
```

## Format

A data frame with 276 rows and 20 variables:

- id:

  case number

- time:

  number of days between registration and the earlier of death,
  transplantion, or study analysis in July, 1986

- status:

  status at endpoint, 0 for censored or transplant, 1 for dead

- trt:

  randomized treatment group: D-penicillmain or placebo

- age:

  in years

- sex:

  m/f

- ascites:

  presence of ascites

- hepato:

  presence of hepatomegaly or enlarged liver

- spiders:

  blood vessel malformations in the skin

- edema:

  0 no edema, 0.5 untreated or successfully treated, 1 edema despite
  diuretic therapy

- bili:

  serum bilirubin (mg/dl)

- chol:

  serum cholesterol (mg/dl)

- albumin:

  serum albumin (g/dl)

- copper:

  urine copper (ug/day)

- alk.phos:

  alkaline phosphotase (U/liter)

- ast:

  aspartate aminotransferase, once called SGOT (U/ml)

- trig:

  triglycerides (mg/dl)

- platelet:

  platelet count

- protime:

  standardized blood clotting time

- stage:

  histologic stage of disease (needs biopsy)

## Source

T Therneau and P Grambsch (2000), Modeling Survival Data: Extending the
Cox Model, Springer-Verlag, New York. ISBN: 0-387-98784-3.

## Details

1.  removed rows with missing data

2.  converted `status` into 0 for censor or transplant, 1 for dead

3.  converted `stage` into an ordered factor.

4.  converted `trt`, `ascites`, `hepato`, `spiders`, and `edema` into
    factors.
