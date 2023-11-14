
test_that(
 desc = "linreg_fit with weights approximately equal to glm()",
 code = {

  nrows <- 1000
  ncols <- 20

  X <- matrix(data = rnorm(nrows*ncols), nrow = nrows, ncol = ncols)

  # X <- cbind(1, X)

  colnames(X) <- c(
   # "intercept",
   paste0("x", seq(ncols))
  )

  Y <- matrix(rnorm(nrows), ncol = 1)

  glm_data <- as.data.frame(cbind(y=as.numeric(Y), X))

  # Fit logistic regression using the custom function

  W <- sample(1:3, nrow(X), replace=TRUE)

  cpp = linreg_fit_exported(X, Y, W, do_scale = FALSE,
                            epsilon = 1e-9, iter_max = 20)

  cpp_scale = linreg_fit_exported(X, Y, W, do_scale = TRUE,
                                  epsilon = 1e-9, iter_max = 20)

  R = lm(y ~ ., weights = as.integer(W), data = glm_data)

  R_summary <- summary(R)

  R_beta_est <- as.numeric(R_summary$coefficients[-1, 'Estimate'])
  R_beta_pvalues <- as.numeric(R_summary$coefficients[-1, 'Pr(>|t|)'])

  expect_equal(cpp[,1], R_beta_est, tolerance = 1e-9)
  # scaling changes estimates by a little (that's expected) but not too much
  expect_equal(cpp_scale[,1], R_beta_est, tolerance = .1)

  expect_equal(cpp[,2], R_beta_pvalues, tolerance = 1e-9)
  # scaling does not change p-values at all
  expect_equal(cpp_scale[,2], R_beta_pvalues, tolerance = 1e-9)


 }
)

test_that(
 desc = 'constant columns dont break linreg',
 code = {

  nrows <- 100
  ncols <- 2

  X <- matrix(data = rnorm(nrows*ncols), nrow = nrows, ncol = ncols)

  # X <- cbind(1, X)

  colnames(X) <- c(
   # "intercept",
   paste0("x", seq(ncols))
  )

  Y <- matrix(rnorm(nrows), ncol = 1)

  W <- rep(1, nrow(X))

  X[,1] <- 1

  zeros <- matrix(0, ncol = 2, nrow = ncols)

  results <- linreg_fit_exported(x_node = X, y_node = Y, w_node = W,
                                 do_scale = T, epsilon = 1e-9, iter_max = 20)

  expect_equal(zeros, results)


 }
)
