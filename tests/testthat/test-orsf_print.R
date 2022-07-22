

test_that("standard object prints fine", {

 object <- orsf(pbc_orsf, time+status ~ . - id)
 expect_invisible(p <- capture.output(print(object)))

})
