
library(survival)

y <- flchain_y

set.seed(329)
weights <- sample(1:5, length(y), replace = TRUE)

bcj <- leaf_surv(y, weights)

kap <- survfit(y ~ 1, weights = weights)

test_that(
 desc = 'leaf_surv has same length as survfit',
 code = {expect_equal(length(kap$surv), nrow(bcj))}
)

test_that(
 desc = 'leaf_surv has same values as survfit',
 code = {expect_lt(max(abs(kap$surv) - bcj[,2]), expected = 1e-10)}
)

bcj_small <- leaf_surv_small(y, weights)

yy <- as.data.frame(as.matrix(y))

yy <- unique(subset(yy, status == 1))

test_that(
 desc = 'leaf_surv_small has same length as survfit (kinda)',
 code = {expect_equal(nrow(yy), nrow(bcj_small))}
)

dt_big <- as.data.table(bcj)

dt_small <- as.data.table(bcj_small)

names(dt_big) <- names(dt_small) <- c('time', 'prob')

dt_big_rdcd <- dt_big[time %in% dt_small$time, ]

test_that(
 desc = 'leaf_surv_small has same values as leaf_surv',
 code = {
  expect_lt(max(abs(dt_big_rdcd$prob-dt_small$prob)), expected = 1e-10)
 }
)
