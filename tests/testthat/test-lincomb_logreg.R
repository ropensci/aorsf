
test_that(
 desc = "logreg_fit with weights approximately equal to glm()",
 code = {

  nrows <- 1000
  ncols <- 20

  set.seed(329)
  X <- matrix(data = rnorm(nrows*ncols), nrow = nrows, ncol = ncols)

  colnames(X) <- c(paste0("x", seq(ncols)))

  Y <- matrix(rbinom(nrows, size = 1, prob = 0.7), ncol = 1)

  W <- sample(1:5, nrow(X), replace = TRUE)

  glm_data <- as.data.frame(cbind(y=as.numeric(Y), X))

  # Fit logistic regression using the custom function

  control <- glm.control()


  glm_fit = glm(y ~ .,
                weights = as.integer(W),
                control = control,
                data = glm_data,
                family = 'binomial')

  .summary <- summary(glm_fit)

  R_estimates <- list(
   beta_est = as.numeric(.summary$coefficients[-1, 'Estimate']),
   beta_pvalues = as.numeric(.summary$coefficients[-1, 'Pr(>|z|)'])
  )


  cpp = logreg_fit_exported(X, Y, W, do_scale = FALSE,
                            epsilon = control$epsilon,
                            iter_max = control$maxit)

  cpp_scaled = logreg_fit_exported(X, Y, W, do_scale = TRUE,
                                   epsilon = control$epsilon,
                                   iter_max = control$maxit)

  expect_equal(cpp[,1], R_estimates$beta_est, tolerance = 1e-3)
  expect_equal(cpp[,2], R_estimates$beta_pvalues, tolerance = 1e-3)
  expect_equal(cpp_scaled[,1], cpp[, 1], tolerance = .1)
  expect_equal(cpp_scaled[,2], cpp[, 2], tolerance = 1e-3)


 }
)

test_that(
 desc = "logreg_fit with rank deficient matrix returns all 0",
 code = {

  X <- cbind(
   c(1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1),
   c(1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1)
  )

  X <- cbind(X, X[,1] + X[,2])

  y <- matrix(rbinom(n = nrow(X), size = 1, prob = 1/2), ncol = 1)

  cpp <- logreg_fit_exported(X, y, rep(1, nrow(X)), FALSE, 1e-9, 20)

  expect_true(all(cpp == 0))

 }
)

