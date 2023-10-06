
test_that(
 desc = 'leaf node stats have same time/surv/chaz as survfit',
 code = {

  for(i in seq_along(mat_list_surv)){

   y <- mat_list_surv[[i]]$y
   w <- mat_list_surv[[i]]$w
   r <- sprout_node_survival_exported(y, w)

   aorsf_surv <- r$prob[[1]]
   aorsf_chaz <- r$chaz[[1]]
   aorsf_time <- r$indx[[1]]
   aorsf_mort <- r$mort[[1]]

   kap_fit <- survfit(Surv(time, status) ~ 1,
                      data = as.data.frame(y),
                      weights = w)

   kap_data <- data.frame(time = kap_fit$time,
                          surv = kap_fit$surv,
                          cumhaz = kap_fit$cumhaz,
                          n_event = kap_fit$n.event)

   kap_data <- subset(kap_data, n_event > 0)

   expect_equal(kap_data$time,   aorsf_time, tolerance = 1e-9)
   expect_equal(kap_data$surv,   aorsf_surv, tolerance = 1e-9)
   expect_equal(kap_data$cumhaz, aorsf_chaz, tolerance = 1e-9)
   expect_equal(sum(kap_data$cumhaz), aorsf_mort, tolerance = 1e-9)

  }

 }
)









