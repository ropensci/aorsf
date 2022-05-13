

# handle 'units' variables
suppressMessages(library(units))

pbc_units <- pbc_orsf

units(pbc_units$time) <- 'days'
units(pbc_units$age) <- 'years'
units(pbc_units$bili) <- 'mg/dl'

test_that("output has expected items", {

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

pbc_units_badclass <- pbc_units
class(attr(pbc_units_badclass$bili, 'units')) <- 'bad_units'


test_that('only symbolic units are allowed', {

 expect_error(unit_info(pbc_units_badclass, 'bili'), 'symbolic_units')

})




