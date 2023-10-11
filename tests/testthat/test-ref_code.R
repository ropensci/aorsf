
fi <- fctr_info(data = pbc_orsf, .names = names(pbc_orsf))

pbc_refcoded <- ref_code(pbc_orsf,
                         fi,
                         names_x_data = c("age", "sex", "stage"))

test_that(
 desc = "reference coding names and types",
 code = {
  expect_named(pbc_refcoded, expected = c("age", "sex_f", "stage"))
  expect_type(pbc_refcoded$stage, 'integer')
  expect_type(pbc_refcoded$sex_f, 'integer')
  expect_type(pbc_refcoded$age, 'double')
 }
)
