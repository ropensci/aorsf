

gini_impurity <- function (vals) {
 counts <- table(vals)
 total <- sum(counts)
 return(sum((counts/total) * (1 - counts/total)))
}

test_that(
 desc = "gini index matches expected answer",
 code = {

  n <- 100

  y <- matrix(rbinom(n, size = 1, prob = 1/2), ncol = 1)
  w <- rep(1, n)
  g <- rbinom(n, size = 1, prob = 1/3)
  vals_1 = factor(y[g==1], levels = c(0,1), labels = c("blue", "green"))
  vals_0 = factor(y[g==0], levels = c(0,1), labels = c("blue", "green"))

  gini_1 <- gini_impurity(vals = vals_1)
  gini_0 <- gini_impurity(vals = vals_0)

  target <- gini_1 * sum(w[g==1]) / sum(w) + gini_0 * sum(w[g==0]) / sum(w)

  cpp <- compute_gini_exported(y, w, g)

  expect_equal(target, cpp)

 }
)

