


test_that(
 desc = 'is_empty true for empty, false for non-empty',
 code = {
  expect_true(is_empty(c()))
  expect_false(is_empty(c(1)))
  expect_false(is_empty(NA_real_))
 }
)

test_that(
 desc = 'is_error true for error, false for non-error',
 code = {
  expect_true(is_error(try("a" + 1, silent = TRUE)))
  expect_false(is_error(1))
 }
)

test_that(
 desc = 'list initiation matches base R list creation',
 code = {
  test_list <- vector(mode = 'list', length = 3)
  names(test_list) <- c("a", "b", "c")
  expect_equal(test_list, list_init(.names = c("a", "b", "c")))
 }
)

test_that(desc = 'last value is last value',
          code = {expect_true(last_value(c(1,2,3)) == 3)})

test_that(desc = 'inheritance for aorsf',
          code = {expect_false(is_aorsf('a'))})





