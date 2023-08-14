
#' @srrstats {G5.4} **Correctness tests** *test that statistical algorithms produce expected results to some fixed test data sets. I use the pbc data and compare the aorsf newton raphson algorithm to the same algorithm used in survival::coxph().*

#' @srrstats {G5.4b} *Correctness tests include tests against previous implementations, explicitly calling those implementations in testing.*

#' @srrstats {G5.0} *tests use the pbc data and the flchain data, two standard datasets in the survival package that are widely studied. The pbc data are also featured in another R package for random forests, i.e., randomForestSRC*


#' @srrstats {G5.6} **Parameter recovery tests** *The coxph newton-raphson algorithm returns coefficient values. The aorsf version matches those within a specified tolerance*

#' @srrstats {ML7.7} *explicitly test optimization algorithms for accuracy. I do not test multiple optimization algorithms because this is the only one I have programmed in aorsf. The optimization algorithm used in coxnet has been thoroughly tested by the glmnet developers.*

#' @srrstats {ML7.10} *The successful extraction of information on paths taken by optimizers is tested.*

iter_max = 20
control <- survival::coxph.control(iter.max = iter_max, eps = 1e-8)

run_cph_test <- function(x, y, method){

 wts <- sample(seq(1:2), size = nrow(x), replace = TRUE)

 tt = survival::coxph.fit(x = x,
                          y = y,
                          strata = NULL,
                          offset = NULL,
                          init = rep(0, ncol(x)),
                          control = control,
                          weights = wts,
                          method = if(method == 0) 'breslow' else 'efron',
                          rownames = NULL,
                          resid = FALSE,
                          nocenter = c(0))

 tt_fit <- survival::coxph(y~x,
                           weights = wts,
                           control = control,
                           ties = if(method == 0) 'breslow' else 'efron')

 tt_inf <- summary(tt_fit)$coefficients[,'Pr(>|z|)']

 xx <- x[, , drop = FALSE]

 bcj = coxph_fit_exported(xx,
                          y,
                          wts,
                          method = method,
                          cph_eps = 1e-8,
                          cph_iter_max = iter_max)

 beta <- bcj$beta

 rownames(beta) <- names(tt$coefficients)
 beta_vec <- beta[, 1, drop = TRUE]

 perc_diff <- function(a,b) abs(a-b) / (abs(0.001 + a+b)/2)

 # maximum percent difference
 max(perc_diff(tt$coefficients, beta_vec))

}


# pbc data ----------------------------------------------------------------

.pbc <- pbc_orsf[order(pbc_orsf$time), ]

.pbc$trt <- as.numeric(.pbc$trt)
.pbc$ascites <- as.numeric(.pbc$ascites)
.pbc$hepato <- as.numeric(.pbc$hepato)

x <- as.matrix(.pbc[, c('trt','age','ascites','hepato','bili')])
y <- survival::Surv(.pbc$time, .pbc$status)

#' @srrstats {G5.6a} *succeed within a defined tolerance rather than recovering exact values.*

test_that(
 desc = 'similar answers for pbc data',
 code = {
  expect_true( run_cph_test(x, y, method = 0) < 1e-2 )
  expect_true( run_cph_test(x, y, method = 1) < 1e-2 )
 }
)

# flchain data ------------------------------------------------------------

data("flchain", package = 'survival')

df <- na.omit(flchain)

df$chapter <- NULL

time <- 'futime'
status <- 'death'

df_nomiss <- na.omit(df)

df_sorted <- df_nomiss[order(df_nomiss[[time]]),]

df_x <- df_sorted
df_x[[time]] <- NULL
df_x[[status]] <- NULL

flchain_x <- model.matrix(~.-1, data = df_x)

flchain_y <- survival::Surv(time = df_sorted[[time]],
                  event = df_sorted[[status]])

x <- flchain_x[, c('age', 'sexF','sample.yr', 'kappa', 'lambda')]
y <- flchain_y

#' @srrstats {G5.6a} *succeed within a defined tolerance rather than recovering exact values.*

test_that(
 desc = 'similar answers for flchain data',
 code = {
  expect_true( run_cph_test(x, y, method = 0) < 1e-2 )
  expect_true( run_cph_test(x, y, method = 1) < 1e-2 )
 }
)



# speed comparison --------------------------------------------------------

data("flchain", package = 'survival')

df <- na.omit(flchain)

df$chapter <- NULL

time <- 'futime'
status <- 'death'

df_nomiss <- na.omit(df)

df_sorted <- df_nomiss[order(df_nomiss[[time]]),]

df_x <- df_sorted
df_x[[time]] <- NULL
df_x[[status]] <- NULL

flchain_x <- model.matrix(~.-1, data = df_x)

flchain_y <- survival::Surv(time = df_sorted[[time]],
                            event = df_sorted[[status]])

x <- flchain_x[, c('age', 'sexF','sample.yr', 'kappa', 'lambda')]
y <- flchain_y

wts <- sample(seq(1:2), size = nrow(x), replace = TRUE)

method = 0

control <- survival::coxph.control(iter.max = 1, eps = 1e-8)

microbenchmark::microbenchmark(

 tt = survival::coxph.fit(x = x,
                          y = y,
                          strata = NULL,
                          offset = NULL,
                          init = rep(0, ncol(x)),
                          control = control,
                          weights = wts,
                          method = if(method == 0) 'breslow' else 'efron',
                          rownames = NULL,
                          resid = FALSE,
                          nocenter = c(0)),

 bcj = coxph_fit_exported(x[, , drop = FALSE],
                          y,
                          wts,
                          method = method,
                          cph_eps = 1e-8,
                          cph_iter_max = control$iter.max)

)

