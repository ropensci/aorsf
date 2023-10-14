

x_rows <- sample(nrow(pbc_mats$x), size = 100)
x_cols <- sample(ncol(pbc_mats$x), size = 10)

beta <- runif(n = length(x_cols))

test_that(
 desc = "submatrix multiplication is correct",
 code = {

  data_cpp_answer <- x_submat_mult_beta_exported(x = pbc_mats$x,
                                                 y = pbc_mats$y,
                                                 w = pbc_mats$w,
                                                 x_rows = x_rows - 1,
                                                 x_cols = x_cols - 1,
                                                 beta = beta)

  target <- pbc_mats$x[x_rows, x_cols] %*% beta

  expect_equal(data_cpp_answer, target)

 })

