
library(survival)

# method = 0 for breslow, 1 for efron

set.seed(32987) # random tests could break by chance

run_cph_test <- function(x, y, method, pval_max = 1/2){

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

 xx <- x[, , drop = FALSE]

 bcj = newtraph_cph_testthat(xx,
                             y,
                             wts,
                             method = method,
                             eps = 1e-8,
                             iter_max = 20,
                             pval_max = pval_max)

 rownames(bcj) <- names(tt$coefficients)
 bcj_vec <- bcj[, 1, drop = TRUE]


 perc_diff <- function(a,b) abs(a-b) / (abs(0.001 + a+b)/2)

 tt$coefficients[tt_inf > pval_max] <- 0

 # maximum percent difference
 max(perc_diff(tt$coefficients, bcj_vec))

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
  expect_true( run_cph_test(x, y, method = 0) < 1e-5 )
  expect_true( run_cph_test(x, y, method = 1) < 1e-5 )
 }
)

# flchain data ------------------------------------------------------------

x <- flchain_x[, c('age', 'sexF','sample.yr', 'kappa', 'lambda')]
y <- flchain_y

test_that(
 desc = 'similar answers for flchain data',
 code = {
  expect_true( run_cph_test(x, y, method = 0) < 1e-5 )
  expect_true( run_cph_test(x, y, method = 1) < 1e-5 )
 }
)
