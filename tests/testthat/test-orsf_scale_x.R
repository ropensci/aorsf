
x_mat <- as.matrix(pbc_orsf[, c('bili',
                                'age',
                                'protime',
                                'chol')])

head(x_mat)

x_scaled <- orsf_scale_cph(x_mat)

test_that(
 desc = 'dimensions and names are conserved',
 code = {
  expect_equal(dim(x_mat), dim(x_scaled))
  expect_equal(colnames(x_mat), colnames(x_scaled))
 }
)


test_that(
 desc = 'transforms are stored as an attribute',
 code = {
  expect_false(is.null(attr(x_scaled, 'transforms')))
  expect_type(attr(x_scaled, 'transforms'), 'double')
 }
)

x_unscaled <- orsf_unscale_cph(x_scaled)

test_that(
 desc = 'inverse transform is equivalent to initial data',
 code = {
  expect_lt(max(abs(x_mat - x_unscaled)), 1e-8)
 }
)

bad_mat <- as.matrix(pbc_orsf[, c('sex',
                                  'age',
                                  'protime',
                                  'chol')])

bad_vec_1 <- c(1)
bad_vec_2 <- rep(-1, nrow(x_mat))
bad_vec_3 <- rep('a', nrow(x_mat))

test_that(
 desc = 'unexpected inputs trigger informative error',
 code = {

  expect_error(orsf_scale_cph(bad_mat), regexp = 'should have type')
  expect_error(orsf_scale_cph(pbc_orsf),'should inherit from')
  expect_error(orsf_scale_cph(x_mat, bad_vec_1), 'length')
  expect_error(orsf_scale_cph(x_mat, bad_vec_2), 'should be > 0')
  expect_error(orsf_scale_cph(x_mat, bad_vec_3), 'should have type')
  expect_error(orsf_unscale_cph(x_mat), 'transforms')

 }
)
