
# R version written using matrixStats

weightedVar <- function (x, w = NULL, idxs = NULL, na.rm = FALSE, center = NULL,
                         ...) {
 n <- length(x)
 if (is.null(w)) {
  w <- rep(1, times = n)
 }
 else if (length(w) != n) {
  stop(sprintf("The number of elements in arguments '%s' and '%s' does not match: %.0f != %.0f",
               "w", "x", length(w), n))
 }
 else if (!is.null(idxs)) {
  w <- w[idxs]
 }
 if (!is.null(idxs)) {
  x <- x[idxs]
  n <- length(x)
 }
 na_value <- NA
 storage.mode(na_value) <- storage.mode(x)
 tmp <- (is.na(w) | w > 0)
 if (!all(tmp)) {
  x <- .subset(x, tmp)
  w <- .subset(w, tmp)
  n <- length(x)
 }
 tmp <- NULL
 if (na.rm) {
  keep <- which(!is.na(x))
  x <- .subset(x, keep)
  w <- .subset(w, keep)
  n <- length(x)
  keep <- NULL
 }

 tmp <- is.infinite(w)
 if (any(tmp)) {
  keep <- tmp
  x <- .subset(x, keep)
  n <- length(x)
  w <- rep(1, times = n)
  keep <- NULL
 }
 tmp <- NULL
 if (n <= 1L)
  return(na_value)
 wsum <- sum(w)
 if (is.null(center)) {
  center <- sum(w * x)/wsum
 }
 x <- x - center
 x <- x^2
 lambda <- 1/(wsum - 1)
 sigma2 <- lambda * sum(w * x)
 x <- w <- NULL

 sigma2
}

var_reduction_R <- function(y, w, g){
  (sum(w) - 1)/sum(w) * weightedVar(y, w = w) -
    (sum(w*g) - 1)/(sum(w))*weightedVar(y, w = w, idxs = which(g == 1)) -
    (sum(w*(1-g)) - 1)/(sum(w))*weightedVar(y, w = w, idxs = which(g == 0))
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
