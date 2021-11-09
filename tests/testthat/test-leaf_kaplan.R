
library(survival)
devtools::load_all()

y <- flchain_y
ymat <- as.matrix(y)

set.seed(329)
weights <- sample(1:5, length(y), replace = TRUE)

bcj <- leaf_kaplan_testthat(ymat, weights)

kap <- survfit(y ~ 1, weights = weights)

kap <- data.frame(n.event = kap$n.event,
                  time = kap$time,
                  surv = kap$surv)

kap <- subset(kap, n.event > 0)

test_that(
 desc = 'leaf_kaplan has same length as survfit',
 code = {expect_equal(nrow(kap), nrow(bcj))}
)

test_that(
 desc = 'leaf_kaplan has same time values as survfit',
 code = {expect_lt(max(abs(kap$time) - bcj[,1]), expected = 1e-10)}
)

test_that(
 desc = 'leaf_kaplan has same surv values as survfit',
 code = {expect_lt(max(abs(kap$surv) - bcj[,2]), expected = 1e-10)}
)



