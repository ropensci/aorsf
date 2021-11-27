
iter_max = 20
control <- survival::coxph.control(iter.max = iter_max, eps = 1e-8)

run_cph_test <- function(x, y, method, pval_max = 1){

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

 bcj = newtraph_cph_testthat(xx,
                             y,
                             wts,
                             method = method,
                             cph_eps_ = 1e-8,
                             iter_max = iter_max,
                             pval_max = pval_max)

 rownames(bcj) <- names(tt$coefficients)
 bcj_vec <- bcj[, 1, drop = TRUE]


 perc_diff <- function(a,b) abs(a-b) / (abs(0.001 + a+b)/2)

 tt$coefficients[tt_inf > pval_max] <- 0

 # maximum percent difference
 max(perc_diff(tt$coefficients, bcj_vec))

}

# pbc data ----------------------------------------------------------------

.pbc <- pbc_orsf[order(pbc_orsf$time), ]

x <- as.matrix(.pbc[, c('trt','age','ascites','hepato','bili')])
y <- survival::Surv(.pbc$time, .pbc$status)

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

test_that(
 desc = 'similar answers for flchain data',
 code = {
  expect_true( run_cph_test(x, y, method = 0) < 1e-2 )
  expect_true( run_cph_test(x, y, method = 1) < 1e-2 )
 }
)
