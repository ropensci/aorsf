
library(survival)

# method = 0 for breslow, 1 for efron

set.seed(32987) # random tests could break by chance

run_cph_test <- function(x, y, method){

 wts <- sample(seq(1:2), size = nrow(x), replace = TRUE)

 tt = coxph.fit(x = x,
                y = y,
                strata = NULL,
                offset = NULL,
                init = rep(0, ncol(x)),
                control = coxph.control(iter.max = 20, eps = 1e-8),
                weights = wts,
                method = if(method == 0) 'breslow' else 'efron',
                rownames = NULL,
                resid = FALSE,
                nocenter = c(0))

 tt_inf <- summary(
   coxph(y~x,
         weights = wts,
         ties = if(method == 0) 'breslow' else 'efron')
   )$coefficients[,'Pr(>|z|)']

 xx <- x[]

 x_transforms = x_scale_wtd(xx, wts)

 bcj = newtraph_cph(xx,
                    y,
                    wts,
                    x_transforms,
                    method = method,
                    eps = 1e-8,
                    iter_max = 20,
                    rescale = TRUE)

 perc_diff <- function(a,b) abs(a-b) / (abs(a+b)/2)

 # maximum percent difference
 coef = max(perc_diff(tt$coefficients, bcj[,1]))

 se = max(perc_diff(sqrt(diag(tt$var)), bcj[,2]))

 pv = max(abs((tt_inf-bcj[,3])))

 return(c(coef = coef, se = se, pv = pv))

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
  expect_true( all(run_cph_test(x, y, method = 0) < 1e-5) )
  expect_true( all(run_cph_test(x, y, method = 1) < 1e-5) )
 }
)

# flchain data ------------------------------------------------------------

x <- flchain_x[, c('age', 'sexF','sample.yr', 'kappa', 'lambda')]
y <- flchain_y

test_that(
 desc = 'similar answers for flchain data',
 code = {
  expect_true( all(run_cph_test(x, y, method = 0) < 1e-5) )
  expect_true( all(run_cph_test(x, y, method = 1) < 1e-5) )
 }
)

# benchmark (leave commented out) -----------------------------------------

# .pbc <-  pbc[order(pbc$time), ]
# .pbc <- .pbc[complete.cases(.pbc), ]
# .pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1
#
# x <- as.matrix(.pbc[, c(5), drop = FALSE])
# y <- Surv(.pbc$time, .pbc$status)
#
#
# wts <- sample(c(1,2,3), nrow(x), replace = TRUE)
#
# doscale <- rep(TRUE, ncol(x))
#
# microbenchmark::microbenchmark(
#
#   survival = coxph.fit(x = x,
#                      y = y,
#                      strata = NULL,
#                      offset = NULL,
#                      init = rep(0, ncol(x)),
#                      control = coxph.control(eps = 1e-8,
#                                              iter.max = 3),
#                      weights = wts,
#                      method = 1,
#                      rownames = NULL,
#                      resid = FALSE,
#                      nocenter = c(0)),
#
#
#
#   orsf2 = {
#     xx = x
#     x_transforms = x_scale_wtd(xx, wts, doscale)
#     newtraph_cph(xx, y, wts, x_transforms, 1, 1e-8,
#                  iter_max = 3, rescale = TRUE)
#   }
# )
