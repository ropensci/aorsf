
test_that(
 desc = 'C-statistic (survival) is close to survival::concordance',
 code = {

  y_mat <- as.matrix(pbc_orsf[, c('time', 'status')])

  sorted <-
   collapse::radixorder(y_mat[, 'time'], -y_mat[, 'status'])

  y_mat <- y_mat[sorted, ]

  x_vars <- c("bili", "chol", "trig")

  n_runs <- 100

  diffs_vec <- diffs_uvec<- diffs_vec_wtd <- diffs_uvec_wtd <-
   vector(mode = 'numeric', length = n_runs)

  for(i in seq(n_runs)){

   x_var <- sample(x_vars, 1)

   s_vec <- pbc_orsf[[x_var]][sorted] + rnorm(nrow(y_mat), mean=0, sd=1)
   s_uvec <- as.integer(s_vec > 1)

   w_1s <- rep(1, nrow(y_mat))
   w_rs <- sample(1:5, nrow(y_mat), replace = TRUE)

   survival_vec     <- 1-oobag_c_survival(y_mat, w_1s, s_vec)
   survival_vec_wtd <- 1-oobag_c_survival(y_mat, w_rs, s_vec)

   survival_uvec     <- 1-oobag_c_survival(y_mat, w_1s, s_uvec)
   survival_uvec_wtd <- 1-oobag_c_survival(y_mat, w_rs, s_uvec)

   aorsf_vec <- compute_cstat_surv_exported_vec(y_mat, w_1s, s_vec, TRUE)
   aorsf_vec_wtd <- compute_cstat_surv_exported_vec(y_mat, w_rs, s_vec, TRUE)

   aorsf_uvec <- compute_cstat_surv_exported_uvec(y_mat, w_1s, s_uvec, TRUE)
   aorsf_uvec_wtd <- compute_cstat_surv_exported_uvec(y_mat, w_rs, s_uvec, TRUE)

   diffs_vec[i] <- abs(survival_vec - aorsf_vec)
   diffs_uvec[i] <- abs(survival_uvec - aorsf_uvec)

   diffs_vec_wtd[i] <- abs(survival_vec_wtd - aorsf_vec_wtd)
   diffs_uvec_wtd[i] <- abs(survival_uvec_wtd - aorsf_uvec_wtd)

  }

  # unweighted is basically identical to cstat from survival
  expect_lt(mean(diffs_vec), 0.001)
  expect_lt(mean(diffs_uvec), 0.001)
  # weighted is very close
  expect_lt(mean(diffs_vec_wtd), 0.01)
  expect_lt(mean(diffs_uvec_wtd), 0.01)

 }
)

test_that(
 desc = 'C-statistic (classification) is close to Hmisc::somers2',
 code = {

  n <- 1000

  p <- list(
   catg = round(rnorm(n), 1),
   ctns = rnorm(n),
   bnry = rbinom(n, size = 1, prob = 1/2)
  )

  y <- matrix(rbinom(n, size = 1, prob = 1/2), ncol = 1)
  w <- sample(1:5, n, replace = TRUE)

  for(i in seq_along(p)){

   target <- Hmisc::somers2(p[[i]], y = y, weights = w)['C']
   answer <- compute_cstat_clsf_exported(y, w, p[[i]])

   expect_equal(as.numeric(target), answer, tolerance = 1e-7)

  }

 }
)
