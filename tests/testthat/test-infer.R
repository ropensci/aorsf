
test_that(
 desc = 'inferred pred horizon is correct',
 code = {

  expect_equal(
   infer_pred_horizon(fit_standard_pbc$fast,
                      pred_type = 'risk',
                      pred_horizon = NULL),
   get_oobag_pred_horizon(fit_standard_pbc$fast)
  )

  expect_equal(
   infer_pred_horizon(fit_standard_pbc$fast,
                      pred_type = 'risk',
                      pred_horizon = 100),
   100
  )

  expect_equal(
   infer_pred_horizon(fit_standard_pbc$fast,
                      pred_type = 'mort',
                      pred_horizon = 100),
   1
  )

  fit_renegade <- fit_standard_pbc$fast

  fit_renegade$pred_horizon <- NULL
  attr(fit_renegade, 'oobag_pred_horizon') <- NULL

  expect_error(infer_pred_horizon(fit_renegade,
                                  pred_type = 'risk',
                                  pred_horizon = NULL),
               regexp = 'could not be found')

 }
)


test_that(
 desc = 'inferred outcome type is correct',
 code = {

  pbc$surv_y <- Surv(pbc_orsf$time, pbc_orsf$status)

  expect_equal(infer_outcome_type(c('time', 'status'), pbc), 'survival')
  expect_equal(infer_outcome_type('surv_y', pbc), 'survival')
  expect_equal(infer_outcome_type('age', pbc), 'regression')
  expect_equal(infer_outcome_type('sex', pbc), 'classification')

 }
)
