
#' @srrstats {G5.4} **Correctness tests** *test that statistical algorithms produce expected results to some fixed test data sets. I use the flchain data and compare the aorsf kaplan meier routine to that of the survival package.*

#' @srrstats {G5.4b} *Correctness tests include tests against previous implementations, explicitly calling those implementations in testing.*

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
#' @srrstats {G5.5} *Correctness tests are run with a fixed random seed*
set.seed(329)

weights <- sample(1:5, length(y), replace = TRUE)

rows <- sort(sample(nrow(ymat), 20))

bcj <- leaf_kaplan_testthat(ymat[rows, ], weights[rows])



kap <- survival::survfit(survival::Surv(ymat[rows,1], ymat[rows,2]) ~ 1,
                         weights = weights[rows])

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
 code = {expect_equal(kap$time, bcj[,1])}
)

test_that(
 desc = 'leaf_kaplan has same surv values as survfit',
 code = {expect_equal(kap$surv, bcj[,2])}
)



