

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

  # won't be used
  pd_x_vals <- matrix(1)
  pd_x_cols <- ncol(pbc_mats$x)

  data_cpp_pd_answer <- x_submat_mult_beta_pd_exported(x = pbc_mats$x,
                                                      y = pbc_mats$y,
                                                      w = pbc_mats$w,
                                                      x_rows = x_rows - 1,
                                                      x_cols = x_cols - 1,
                                                      beta = beta,
                                                      pd_x_vals = pd_x_vals,
                                                      pd_x_cols = pd_x_cols)


  target <- pbc_mats$x[x_rows, x_cols] %*% beta

  expect_equal(data_cpp_answer, target)
  expect_equal(data_cpp_pd_answer, target)

 })

test_that(
 desc = "submatrix multiplication with PD values is correct",
 code = {

  pd_x_vals <- c(0, 0)
  pd_x_cols <- x_cols[1:2] - 1

  data_cpp_pd_answer <- x_submat_mult_beta_pd_exported(x = pbc_mats$x,
                                                      y = pbc_mats$y,
                                                      w = pbc_mats$w,
                                                      x_rows = x_rows - 1,
                                                      x_cols = x_cols - 1,
                                                      beta = beta,
                                                      pd_x_vals = pd_x_vals,
                                                      pd_x_cols = pd_x_cols)

  x_pd <- pbc_mats$x
  x_pd[, x_cols[c(1,2)]] <- 0

  target <- x_pd[x_rows, x_cols] %*% beta

  expect_equal(data_cpp_pd_answer, target)


 }
)

