
y <- as.numeric(penguins_orsf$species)-1
y_expand <- expand_y_clsf_exported(y, n_class = 3)

test_that(
 desc = "unweighted probabilities are correct",
 code = {
  w <- sample(1:5, length(y), replace = TRUE)
  y_probs_wtd <- compute_pred_prob_exported(y_expand, w)
  target_probs_wtd <- apply(y_expand, 2, weighted.mean, w)
  target_probs_wtd <- c(1-sum(target_probs_wtd), target_probs_wtd)
  expect_equal(y_probs_wtd, matrix(target_probs_wtd, ncol = 1))
 }
)

test_that(
 desc = "unweighted probabilities are correct",
 code = {
  w_1 <- rep(1, length(y))
  y_probs_raw <- compute_pred_prob_exported(y_expand, w_1)
  target_probs_raw <- matrix(prop.table(table(y)), ncol = 1)
  expect_equal(y_probs_raw, target_probs_raw)
 }
)

