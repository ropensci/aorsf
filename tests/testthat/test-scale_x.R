weighted_sd <- function(x, w){
 mu <- weighted.mean(x, w)
 w_sum <- sum(w)
 m <- length(w)
 sqrt( sum(w * (x-mu)^2) / ( (m-1) * w_sum / m ) )
}

nrows <- 1000
ncols <- 20

X <- matrix(data = rnorm(nrows*ncols), nrow = nrows, ncol = ncols)

colnames(X) <- paste0("x", seq(ncols))

weights <- sample(1:3, nrow(X), replace = TRUE)
x_means <- apply(X, 2, weighted.mean, w = weights)
x_sds <- apply(X, 2, weighted_sd, w = weights)

res <- scale_x_exported(X, w = weights)

expect_equal(res$transforms[, 1], as.numeric(x_means), tolerance = 1e-9)
expect_equal(res$transforms[, 2], as.numeric(x_sds), tolerance = 1e-9)
