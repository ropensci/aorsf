
library(survival)

data("flchain", package = 'survival')

df <- na.omit(flchain)

df$chapter <- NULL

time <- 'futime'
status <- 'death'

df_nomiss <- na.omit(df)

df_sorted <- df_nomiss[order(df_nomiss[[time]]),]

df_x <- df_sorted
df_x[[time]] <- NULL
df_x[[status]] <- NULL

flchain_x <- model.matrix(~.-1, data = df_x)

flchain_y <- Surv(time = df_sorted[[time]],
                  event = df_sorted[[status]])

y <- flchain_y
ymat <- as.matrix(y)

set.seed(329)
weights <- sample(1:5, length(y), replace = TRUE)

bcj <- leaf_kaplan_testthat(ymat, weights)

# run this code locally to creat kap below
kap <- survival::survfit(survival::Surv(ymat[,1], ymat[,2]) ~ 1,
                         weights = weights)

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



