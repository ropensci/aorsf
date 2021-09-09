
library(survival)

y <- flchain_y

set.seed(329)

weights <- sample(1:5, length(y), replace = TRUE)

g <- rbinom(nrow(y), size = 1, prob = 1/2)

log_rank_test(y, g)

log_rank_test_wtd(y, g, rep(1, length(g)))

survdiff(y~g)$chisq

data_wtd <- data.frame(cbind(y), g=g, weights)

data_wtd_unnest <- data_wtd[rep(seq(nrow(data_wtd)), data_wtd$weights), ]

log_rank_test_wtd(y, g, weights)

survdiff(Surv(time, status) ~ g, data = data_wtd_unnest)$chisq

