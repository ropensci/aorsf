
run_cph_test <- function(x, y, w, method){

 control <- coxph.control(iter.max = 20, eps = 1e-8)

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

 fit <- coxph(y ~ x,
              weights = w,
              control = control,
              method = if(method == 0) 'breslow' else 'efron')

 fit_stats <- as.data.frame(
  summary(fit)$coefficients[,c("coef", "Pr(>|z|)")]
 )

 fit_coefs <- fit_stats$coef
 fit_pvalues <- fit_stats$`Pr(>|z|)`

 xx <- x[, , drop = FALSE]

 bcj = coxph_fit_exported(xx,
                          y,
                          w,
                          method = method,
                          epsilon = control$eps,
                          iter_max = control$iter.max)

 expect_equal(as.numeric(fit_coefs), bcj$beta, tolerance = 1e-5)
 expect_equal(as.numeric(fit_pvalues),  bcj$pvalues, tolerance = 1e-5)

}

for(i in seq_along(mat_list_surv)){

 x <- mat_list_surv[[i]]$x
 y <- mat_list_surv[[i]]$y
 w <- mat_list_surv[[i]]$w

 run_cph_test(x, Surv(y), w, method = 0)
 run_cph_test(x, Surv(y), w, method = 1)

}
