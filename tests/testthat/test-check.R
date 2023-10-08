

test_that(
 desc = 'errors trigger on cue',
 code = {

  expect_error(
   check_arg_uni(uni = c('a', 'b', 'd'),
                 arg_name = 'ok',
                 expected_uni = c('a', 'b', 'c')),
   regexp = 'should contain'
  )

  expect_error(
   check_arg_is_integer(arg_name = 'x',
                        arg_value = c(1, 2, 3.3)),
   regexp = 'only integer values'
  )

  expect_error(
   check_arg_is_integer(arg_name = 'x',
                        arg_value = c(1.3)),
   regexp = 'be an integer value'
  )

  expect_error(
   check_arg_length(arg_name = 'x',
                    arg_value = c(1,2,3),
                    expected_length = 4),
   regexp = 'should have length <4>'
  )

  abc <- data.frame(a=1, b=2, c=3)
  cde <- data.frame(c=3, d=4, e=5)

  expect_error(
   check_new_data_names(new_data = cde,
                        ref_names = names(abc),
                        label_new = 'new data',
                        label_ref = 'training data',
                        check_new_in_ref = TRUE,
                        check_ref_in_new = FALSE),
   regexp = 'd and e'
  )

  expect_error(
   check_new_data_names(new_data = cde,
                        ref_names = names(abc),
                        label_new = 'new data',
                        label_ref = 'training data',
                        check_new_in_ref = FALSE,
                        check_ref_in_new = TRUE),
   regexp = 'a and b'
  )

  expect_error(
   check_new_data_names(new_data = cde,
                        ref_names = names(abc),
                        label_new = 'new data',
                        label_ref = 'training data',
                        check_new_in_ref = TRUE,
                        check_ref_in_new = TRUE),
   regexp = 'Also,'
  )

 }

)



