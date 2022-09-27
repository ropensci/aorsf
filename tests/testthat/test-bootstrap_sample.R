


# test setup --------------------------------------------------------------

# make an r function that matches the function in cpp
r_fun <- function(x, y, wts){

 s = seq(0, 10)
 n_rows = nrow(x)
 probs = dbinom(s, n_rows, 1.0/n_rows, FALSE)

 boot_wts = sample(s, n_rows, TRUE, probs);

 if(length(wts) > 0){
  boot_wts = boot_wts * wts;
 }

 return(matrix(boot_wts, ncol = 1));

}

n <- 1e6
x <- y <- matrix(rnorm(n = n), ncol = 1)

weights_2s = rep(2, nrow(x))
weights_empty = double(0)


# if any object's memory is copied, it will make the test output messy
# (I am not sure how to formally make a test fail when memory is copied)
tracemem(x)
tracemem(y)
tracemem(weights_2s)

set.seed(329)
samp_2s_wts <- bootstrap_sample_testthat(x, y, weights_2s)
set.seed(329)
samp_no_wts <- bootstrap_sample_testthat(x, y, weights_empty)
set.seed(329)
samp_2s_wts_r <- r_fun(x, y, weights_2s)
set.seed(329)
samp_no_wts_r <- r_fun(x, y, weights_empty)


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

test_that(
 desc = "bootstrap weights match (R vs cpp)",
 code = {
  expect_equal(samp_2s_wts, samp_2s_wts_r)
  expect_equal(samp_no_wts, samp_no_wts_r)
 }
)


# test performance --------------------------------------------------------

# TODO: figure out why armadillo seems to run slower than R for this??

# microbenchmark::microbenchmark(
#  r_fun = r_fun(x, y, wts = weights_empty),
#  c_fun = bootstrap_sample_testthat(x, y, weights_empty)
# )
#
# microbenchmark::microbenchmark(
#  r_fun = r_fun(x, y, wts = weights_2s),
#  c_fun = bootstrap_sample_testthat(x, y, weights_2s)
# )


