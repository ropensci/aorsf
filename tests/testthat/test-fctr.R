
pbc_char <- pbc_orsf

pbc_char$sex <- as.character(pbc_char$sex)
pbc_char$trt <- as.character(pbc_char$trt)

test_that("no characters allowed", {

 expect_error(fctr_check(pbc_char, .names = names(pbc_char)),
              'trt and sex')


})
