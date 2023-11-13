


var_reduction_R <- function(y, w, g){
  (sum(w) - 1)/sum(w) * weighted_variance(y, w = w) -
    (sum(w*g) - 1)/(sum(w))*weighted_variance(y, w = w, idxs = which(g == 1)) -
    (sum(w*(1-g)) - 1)/(sum(w))*weighted_variance(y, w = w, idxs = which(g == 0))
}

test_that(
  desc = 'computed variance reduction close to matrixStats::weighted_variance',
  code = {

    n_runs <- 100

    diffs_vec <- vector(mode = 'numeric', length = n_runs)

    for(i in seq(n_runs)){

      y <- rnorm(100)
      w <- runif(100, 0, 2)
      g <- rbinom(100, 1, 0.5)
      diffs_vec[i] <- abs(compute_var_reduction_exported(y, w, g) -
                            var_reduction_R(y, w, g))
    }

    # basically identical to R version
    expect_equal(diffs_vec, rep(0, length(diffs_vec)), tolerance = 1e-6)
  }
)


# microbenchmark::microbenchmark(
#   cpp = compute_var_reduction_exported(y, w, g),
#   r = var_reduction_R(y, w, g),
#   times = 10000
# )
