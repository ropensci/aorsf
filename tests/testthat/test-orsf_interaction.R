
object <- orsf(data_train = pbc_orsf,
               formula = Surv(time, status) ~ . - id,
               n_tree = 200)

test_that(
 'min_pairwise_obs is vetted',
 code = {
  expect_error(
   orsf_interaction(object, min_pairwise_obs = 0), regexp = '>= 1'
  )
  expect_error(orsf_interaction(object, min_pairwise_obs = c(2,3)),
               regexp = 'has length <2>')
 }

)

oi <- orsf_interaction(object, min_pairwise_obs = 1)

test_that("all pairs included", {
  expect_equal(nrow(oi), choose(18, 2))
})

