

train <- sample(nrow(pbc_orsf), size = 170)

aorsf = orsf(formula = time + status  ~ . - id,
             data = pbc_orsf[train, ],
             mtry = 5,
             n_split = 10,
             n_tree = 50,
             oobag_pred = TRUE,
             leaf_min_obs = 15)

new_data <- pbc_orsf[-train, ]

test_that(
 desc = 'predictions are bounded',
 code = {
  p <- predict(aorsf, new_data = new_data, times = 1000)
  expect_true(all(p < 1) && all(p > 0))
 })

bad_data <- new_data
bad_data$trt <- factor(bad_data$trt)

test_that(
 desc = 'malicious data types are detected',
 code = {
  expect_error(
   object = predict(aorsf, bad_data),
   regexp = "\\<trt\\>"
  )
 }
)

bad_data <- new_data
