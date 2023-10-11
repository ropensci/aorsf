
run_cph_test <- function(x, y, w, method){

 control <- coxph.control(iter.max = 20, eps = 1e-8)

 start <- Sys.time()

 tt = survival::coxph.fit(x = x,
                          y = y,
                          strata = NULL,
                          offset = NULL,
                          init = rep(0, ncol(x)),
                          control = control,
                          weights = w,
                          method = if(method == 0) 'breslow' else 'efron',
                          rownames = NULL,
                          resid = FALSE,
                          nocenter = c(0))

 stop <- Sys.time()

 tt_time <- stop-start

 xx <- x[, , drop = FALSE]

 start <- Sys.time()

 bcj = coxph_fit_exported(xx,
                          y,
                          w,
                          method = method,
                          cph_eps = control$eps,
                          cph_iter_max = control$iter.max)

 stop <- Sys.time()

 bcj_time <- stop-start

 expect_equal(as.numeric(tt$coefficients), bcj$beta, tolerance = control$eps)

 expect_equal(diag(tt$var),  bcj$var, tolerance = control$eps)

 # list(bcj_time = bcj_time, tt_time = tt_time)

}

for(i in seq_along(mat_list_surv)){

 x <- mat_list_surv[[i]]$x
 y <- mat_list_surv[[i]]$y
 w <- mat_list_surv[[i]]$w

 run_cph_test(x, Surv(y), w, method = 0)
 run_cph_test(x, Surv(y), w, method = 1)

}
