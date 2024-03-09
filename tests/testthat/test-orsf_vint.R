
test_that(
 desc = "orsf_vint orders interactions correctly for requested predictors",
 code = {

  skip_on_cran()

  set.seed(329)

  n <- 100

  data <- data.frame(
   x1 = rnorm(100),
   x2 = rnorm(100),
   x3 = rnorm(100)
  )

  data$y = with(data, expr = x1 + x2 + x3 + 1/2*x1*x2 + x2*x3 + rnorm(n))

  forest <- orsf(data, y ~ ., n_tree = 100)

  vints_1 <- orsf_vint(forest)
  vints_2 <- orsf_vint(forest, predictors = c("x1", "x3"))

  expect_equal(vints_1$interaction[1], "x2..x3")
  expect_equal(vints_1$interaction[2], "x1..x2")
  expect_equal(vints_1$interaction[3], "x1..x3")
  expect_equal(nrow(vints_2), 1)
  expect_equal(vints_2$score,
               vints_1$score[vints_1$interaction == 'x1..x3'])

 }
)

test_that(
 desc = "vint succeeds on categorical forests",
 code = {

  skip_on_cran()

  fit <- orsf(species ~ ., data = penguins_orsf)

  vints <- orsf_vint(fit)

  expect_true(all(levels(penguins_orsf$species) %in% vints$class))

  penguins_bnry <- penguins_orsf
  penguins_bnry$species <- factor(penguins_bnry$species == "Adelie",
                                  levels = c(FALSE, TRUE),
                                  labels = c("Other", "adelie"))

  fit <- orsf(species ~ ., data = penguins_bnry)

  vints <- orsf_vint(fit)

  expect_true(all(levels(penguins_bnry$species) %in% vints$class))

 }
)
