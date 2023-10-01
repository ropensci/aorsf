
#' @srrstats {G5.4} **Correctness tests** *test that statistical algorithms produce expected results to some fixed test data sets. I simulate arbitrary data and compare the aorsf likelihood ratio test to the same algorithm used in survival::survdiff().*

#' @srrstats {G5.4b} *Correctness tests include tests against previous implementations, explicitly calling those implementations in testing.*

#' @srrstats {G5.5} *Correctness tests are run with a fixed random seed*
set.seed(329555)

#' @srrstats {G5.6} **Parameter recovery tests** *the likelihood ratio test returns expected values consistent with the survival implementation for randomly generated data*

#' @srrstats {ML7.8} *Explicitly test my implementation of the likelihood ratio test, used as a loss function when determining the best split for a given node. I do not test other loss functions because this is the only loss function that aorsf implements.*

library(survival)

y <- pbc_orsf[,c('time', 'status')]

sorted <- order(y[, 1], -y[, 2])

y_sort <- as.matrix(y[sorted, ])

g_sort <- as.integer(pbc_orsf$edema == 1)[sorted]

test_that(
 desc = "log-rank test matches survival package",
 code = {

  # with weights
  w <- sample(0:2, nrow(y_sort), replace = TRUE)
  surv_weights <- rep(seq(nrow(y)), times = w)

  y_surv <- y_sort[surv_weights, ]
  g_surv <- g_sort[surv_weights]

  expect_equal(
   # log-rank test from survival
   survival::survdiff(survival::Surv(y_surv[,1], y_surv[,2]) ~ g_surv)$chisq,
   # log-rank test from aorsf
   compute_logrank_exported(y_sort, w, g_sort),
   tolerance = 1e-9
  )

  # without weights
  w <- rep(1, nrow(y_sort))
  surv_weights <- seq(nrow(y_sort))

  expect_equal(
   # log-rank test from survival
   survival::survdiff(survival::Surv(y_sort[,1], y_sort[,2]) ~ g_sort)$chisq,
   # log-rank test from aorsf
   compute_logrank_exported(y_sort, w, g_sort),
   tolerance = 1e-9
  )

 }
)


# # benchmark does not need to be tested every time
#
# bm <- microbenchmark::microbenchmark(
#  surv = survdiff(Surv(y_surv[, 1], y_surv[, 2]) ~ g_surv)$chisq,
#  aorsf = compute_logrank_exported(y_sort, w, g_sort),
#  times = 50
# )
#
# expect_lt(
#  median(bm$time[bm$expr == 'aorsf']),
#  median(bm$time[bm$expr == 'surv'])
# )



