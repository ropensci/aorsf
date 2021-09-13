
library(survival)

# method = 0 for breslow, 1 for efron

run_cph_test <- function(x, y, method){

 wts <- sample(seq(1:4), size = nrow(x), replace = TRUE)

 doscale <- rep(TRUE, ncol(x))

 tt = coxph.fit(x = x,
                    y = y,
                    strata = NULL,
                    offset = NULL,
                    init = rep(0, ncol(x)),
                    control = coxph.control(iter.max = 20, eps = 1e-8),
                    weights = wts,
                    method = method,
                    rownames = NULL,
                    resid = FALSE,
                    nocenter = c(0))

 xx <- x

 x_transforms = x_scale_cph(xx, wts, doscale)

 bcj = newtraph_cph(xx, y, wts, x_transforms, 0, 1e-8,
                    iter_max = 20, rescale = TRUE)

 perc_diff <- function(a,b) abs(a-b) / (abs(a+b)/2)

 # maximum percent difference
 coef = max(perc_diff(tt$coefficients, bcj$beta))

 se = max(perc_diff(sqrt(diag(tt$var)), bcj$se))

 return(c(coef = coef, se = se))

}

# pbc data ----------------------------------------------------------------


.pbc <-  pbc[order(pbc$time), ]
.pbc <- .pbc[complete.cases(.pbc), ]
.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1

x <- as.matrix(.pbc[, -c(1,2,3,6), drop = FALSE])
y <- Surv(.pbc$time, .pbc$status)

test_that(
 desc = 'similar answers for pbc data',
 code = {
  expect_true( all(run_cph_test(x, y, method = 0) < 1e-10) )
  expect_true( all(run_cph_test(x, y, method = 1) < 1e-10) )
 }
)

# flchain data (with weights) ---------------------------------------------

x <- flchain_x[, c('age', 'sexF','sample.yr', 'kappa', 'lambda')]
y <- flchain_y

test_that(
 desc = 'similar answers for flchain data',
 code = {
  expect_true( all(run_cph_test(x, y, method = 0) < 1e-10) )
  expect_true( all(run_cph_test(x, y, method = 1) < 1e-10) )
 }
)
