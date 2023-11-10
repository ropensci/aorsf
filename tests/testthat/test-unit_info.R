

# handle 'units' variables. All of these tests are skipped
# on CRAN because for some reason when I load and use the
# units package it makes valgrind detect possible memory leaks.

test_that(
 'orsf tracks meta data for units class variables',
 code = {

  # units may have memory leaks
  skip_on_cran()

  suppressMessages(library(units))

  pbc_units <- pbc_orsf

  units(pbc_units$time) <- 'days'
  units(pbc_units$age) <- 'years'
  units(pbc_units$bili) <- 'mg/dl'

  fit_units <- orsf(pbc_units, Surv(time, status) ~ . - id, n_tree=1)

  expect_equal(
   fit_units$get_var_unit('time'),
   list( numerator = "d", denominator = character(0), label = "d")
  )

  expect_equal(
   fit_units$get_var_unit('age'),
   list(numerator = "years", denominator = character(0), label = "years")
  )

  expect_equal(
   fit_units$get_var_unit('bili'),
   list(numerator = "mg", denominator = "dl", label = "mg/dl")
  )

 }

)


test_that("output has expected items", {

 skip_on_cran()

 suppressMessages(library(units))

 pbc_units <- pbc_orsf

 units(pbc_units$time) <- 'days'
 units(pbc_units$age) <- 'years'
 units(pbc_units$bili) <- 'mg/dl'
 ui <- unit_info(pbc_units, c('time', 'age', 'bili'))

 expect_equal(
  ui,
  list(
   time = list(
    numerator = "d",
    denominator = character(0),
    label = "d"
   ),
   age = list(
    numerator = "years",
    denominator = character(0),
    label = "years"
   ),
   bili = list(
    numerator = "mg",
    denominator = "dl",
    label = "mg/dl"
   )
  )
 )

 expect_true(is_empty(unit_info(pbc_units, c())))

})


test_that('only symbolic units are allowed', {


 skip_on_cran()

 suppressMessages(library(units))

 pbc_units <- pbc_orsf

 units(pbc_units$bili) <- 'mg/dl'

 class(attr(pbc_units$bili, 'units')) <- 'bad_units'

 expect_error(unit_info(pbc_units, 'bili'), 'symbolic_units')

})


test_that(
 desc = 'inconsistent units are detected',
 code = {

  skip_on_cran()

  suppressMessages(library(units))

  pbc_units <- pbc_orsf
  units(pbc_units$age) <- 'years'

  pbc_test_units <- pbc_orsf
  units(pbc_test_units$age) <- 'days'

  fit <- orsf(formula = time + status  ~ . - id,
              data = pbc_units,
              n_tree = n_tree_test)

  expect_error(predict(fit, new_data = pbc_test_units, pred_horizon = 1000),
               'has unit')

 }
)
