

y <- matrix(
 data = c(1,1,
          2,0,
          3,1,
          4,0,
          5,1,
          6,1),
 byrow = TRUE,
 ncol = 2
)

# column 1 is constant everywhere            -> return value should be 0
# column 2 is constant where events occurred -> return value should be 0
# column 3 is unique everywhere              -> return value should be 1
# column 4 is unique where events occurred   -> return value should be 1

x <- matrix(
 data = c(1, 2, 3, 1,
          1, 0, 4, 1,
          1, 2, 5, 3,
          1, 0, 6, 3,
          1, 2, 7, 4,
          1, 2, 8, 4),
 byrow = TRUE,
 ncol = 4
)

test_that(
 desc = "constant columns are detected in X matrix",
 code = {
  # include all rows
  expect_equal(
   which_cols_valid_exported(y, x, rows_node = seq(0, nrow(x)-1), mtry = 1),
   matrix(data = c(0, 0, 1, 1), ncol = 1)
  )

  # same deal if you just include event rows
  expect_equal(
   which_cols_valid_exported(y, x, rows_node = c(0, 2, 4, 5), mtry = 1),
   matrix(data = c(0, 0, 1, 1), ncol = 1)
  )

  # all cols should be constant if you remove rows with events
  expect_equal(
   which_cols_valid_exported(y, x, rows_node = c(1, 3), mtry = 1),
   matrix(data = c(0, 0, 0, 0), ncol = 1)
  )
 }
)



