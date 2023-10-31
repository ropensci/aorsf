y <- as.numeric(penguins_orsf$species)-1

test_that(
 desc = "y expands to n_class - 1 column matrix (ref coded)",
 code = {

  y_expand <- expand_y_clsf(y, n_class = 3)
  zeros <- which(y == 0)
  ones <- which(y == 1)
  twos <- which(y == 2)

  expect_true(all(y_expand[zeros, ] == 0))
  # ones should be 1 in column 1, 0 o.w.
  expect_true(all(y_expand[ones, 1] == 1))
  expect_true(all(y_expand[ones, 2] == 0))
  # twos should be 1 in column 2, 0 o.w.
  expect_true(all(y_expand[twos, 1] == 0))
  expect_true(all(y_expand[twos, 2] == 1))

 }
)
