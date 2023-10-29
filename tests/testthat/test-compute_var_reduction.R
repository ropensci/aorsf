
# R version written using matrixStats
var_reduction_R <- function(y, w, g){
  (sum(w) - 1)/sum(w) * matrixStats::weightedVar(y, w = w) - 
    (sum(w*g) - 1)/(sum(w))*matrixStats::weightedVar(y, w = w, idxs = which(g == 1)) -
    (sum(w*(1-g)) - 1)/(sum(w))*matrixStats::weightedVar(y, w = w, idxs = which(g == 0))
}

test_that(
  desc = 'computed variance reduction close to matrixStats::weightedVar',
  code = {
    
    n_runs <- 100
    
    diffs_vec <- vector(mode = 'numeric', length = n_runs)
    
    for(i in seq(n_runs)){
      
      y <- rnorm(100)
      w <- runif(100, 0, 2)
      g <- rbinom(100, 1, 0.5)
      diffs_vec[i] <- abs(compute_var_reduction(y, w, g) - 
                            var_reduction_R(y, w, g))
    }
    
    # unweighted is basically identical to cstat from survival
    expect_lt(mean(diffs_vec), 1e-6)
  }
)


# # The cpp implementation is 80+ times faster than the implementation using
# # matrixStats::weightedVar
# microbenchmark::microbenchmark(
#   cpp = compute_var_reduction(y, w, g),
#   r = var_reduction_R(y, w, g),
#   times = 10000
# )
