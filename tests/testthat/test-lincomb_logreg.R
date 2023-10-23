
test_that(
 desc = "logreg_fit with weights approximately equal to glm()",
 code = {

  nrows <- 1000
  ncols <- 20

  X <- matrix(data = rnorm(nrows*ncols), nrow = nrows, ncol = ncols)

  # X <- cbind(1, X)

  colnames(X) <- c(
   # "intercept",
   paste0("x", seq(ncols))
  )

  Y <- matrix(rbinom(nrows, size = 1, prob = 0.7), ncol = 1)

  glm_data <- as.data.frame(cbind(y=as.numeric(Y), X))

  # Fit logistic regression using the custom function

  W <- sample(1:3, nrow(X), replace=TRUE)

  control <- glm.control()

  cpp = logreg_fit_exported(X, Y, W, do_scale = T,
                            epsilon = control$epsilon,
                            iter_max = control$maxit)

  R = glm(y ~ .,
          weights = as.integer(W),
          control = control,
          data = glm_data,
          family = 'binomial')

  R_summary <- summary(R)

  R_beta_est <- as.numeric(R_summary$coefficients[-1, 'Estimate'])
  R_beta_pvalues <- as.numeric(R_summary$coefficients[-1, 'Pr(>|z|)'])

  expect_equal(cpp[,1], R_beta_est, tolerance = 1e-3)
  expect_equal(cpp[,2], R_beta_pvalues, tolerance = 1e-3)


 }
)

