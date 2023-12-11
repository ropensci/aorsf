
set.seed(329)

data <- data.frame(
 x1 = rnorm(500),
 x2 = rnorm(500),
 x3 = rnorm(500)
)

data$y = with(data, expr = x1 + x2 + x3 + 1/2*x1*x2 + x2*x3 + rnorm(500))

forest <- orsf(data, y ~ ., n_tree = 5)

vints_1 <- orsf_vint(forest)
vints_2 <- orsf_vint(forest, predictors = c("x1", "x3"))

test_that(
 desc = "orsf_vint orders interactions correctly",
 code = {

  expect_equal(vints_1$interaction[1], "x2..x3")
  expect_equal(vints_1$interaction[2], "x1..x2")
  expect_equal(vints_1$interaction[3], "x1..x3")

 }
)

test_that(
 desc = "orsf_vint uses only predictors requested",
 code = {

  expect_equal(nrow(vints_2), 1)

 }
)

test_that(
 desc = "interaction score does not depend on unused predictors",
 code = {
  expect_equal(vints_2$score,
               vints_1$score[vints_1$interaction == 'x1..x3'])
 }
)
