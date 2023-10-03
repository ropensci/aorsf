

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
  # column 1 is constant everywhere            -> return value should be FALSE
  # column 2 is constant where events occurred -> return value should be FALSE
  # column 3 is unique everywhere              -> return value should be TRUE
  # column 4 is unique where events occurred   -> return value should be TRUE
  answers <- c(FALSE, FALSE, TRUE, TRUE)

  for(j in seq(ncol(x))){
   expect_equal(
    is_col_splittable_exported(x, y, r = seq(0, nrow(x)-1), j-1),
    answers[j]
   )
  }

  # same deal if you just include event rows
  for(j in seq(ncol(x))){
   expect_equal(
    is_col_splittable_exported(x, y, r = c(0, 2, 4, 5), j-1),
    answers[j]
   )
  }

  # all cols should be not splittable if you remove rows with events
  for(j in seq(ncol(x))){
   expect_equal(
    is_col_splittable_exported(x, y, r = c(1, 3), j-1),
    FALSE
   )
  }
 }
)



