
test_that(
 desc = 'C-statistic is close to survival::concordance',
 code = {

  y_mat <- as.matrix(pbc_orsf[, c('time', 'status')])
  s_vec <- pbc_orsf$bili[sorted]

  sorted <-
   collapse::radixorder(y_mat[, 'time'],  # order this way for risk sets
                        -y_mat[, 'status']) # order this way for oob C-statistic.


  y_mat <- y_mat[sorted, ]
  s_vec <- s_vec[sorted]

  s_uvec <- as.integer(s_vec > 1)

  survival_vec <- survival::concordancefit(
   y = survival::Surv(y_mat),
   x = 1-s_vec
  )$concordance

  survival_uvec <- survival::concordancefit(
   y = survival::Surv(y_mat),
   x = 1-s_uvec
  )$concordance

  w <- rep(1, nrow(y_mat))

  aorsf_vec <- compute_cstat_exported_vec(
   y_mat, w, s_vec, pred_is_risklike = TRUE
  )

  aorsf_uvec <- compute_cstat_exported_uvec(
   y_mat, w, s_uvec, pred_is_risklike = TRUE
  )

  # close enough to cstat from survival
  expect_lt(abs(survival_vec - aorsf_vec), 0.01)
  expect_lt(abs(survival_uvec - aorsf_uvec), 0.01)



 }
)

# # aorsf about 2 times faster
# microbenchmark::microbenchmark(
#  survival = survival::concordancefit(
#   y = survival::Surv(y_mat),
#   x = s_vec
#  ),
#  aorsf = compute_cstat_exported_vec(
#   y_mat, w, s_vec, pred_is_risklike = FALSE
#  )
# )
#
# # aorsf about 2 times faster
# microbenchmark::microbenchmark(
#  survival = survival::concordancefit(
#   y = survival::Surv(y_mat),
#   x = s_uvec
#  ),
#  aorsf = compute_cstat_exported_uvec(
#   y_mat, w, s_uvec, pred_is_risklike = FALSE
#  )
# )
