


# test setup --------------------------------------------------------------

n <- 1e6
x = matrix(rnorm(n = n), ncol = 1)
y_dbl = x
y_int = sample.int(2, size = n, replace = TRUE)

weights_2s = rep(2, nrow(x))
weights_empty = double(0)

set.seed(329)
samp_2s_wts <- bootstrap_sample_testthat(x, y_dbl, y_int, weights_2s)
set.seed(329)
samp_no_wts <- bootstrap_sample_testthat(x, y_dbl, y_int, weights_empty)


# test behavior -----------------------------------------------------------

test_that(
 desc = "max times sampled is < 10",
 code = {
  expect_lt(object = max(samp_no_wts), expected = 10)
 }
)

test_that(
 desc = "bootstrap weights are multiplied by data weights",
 code = {
  # also shows bootstrap weights are replicable by seed set in R
  expect_equal(samp_no_wts, samp_2s_wts / 2)
 }
)

test_that(
 desc = "bootstrap weights include roughly 63.2% of original sample",
 code = {
  expect_equal(mean(samp_no_wts > 0), 0.632, tolerance = 0.01)
 }
)

# test performance --------------------------------------------------------

# r_fun <- function(x, wts){
#
#  s = seq(0, 10)
#  n_rows = nrow(x)
#  probs = dbinom(s, n_rows, 1.0/n_rows, FALSE)
#
#  boot_wts = sample(s, n_rows, TRUE, probs);
#
#  if(length(wts) > 0){
#   boot_wts = boot_wts * wts;
#  }
#  return(boot_wts);
#
# }
#
# microbenchmark::microbenchmark(
#  r_fun = r_fun(x, wts = weights_empty),
#  c_fun = bootstrap_sample_testthat(x, y_dbl, y_int, weights_empty)
# )


