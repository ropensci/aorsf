

test_that("no errors for standard object", {

 object <- orsf(pbc_orsf, time+status ~ . - id)
 expect_invisible(print(object))

})
