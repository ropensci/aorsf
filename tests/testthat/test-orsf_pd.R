test_that(
 desc = "oob stops if there are no data",
 code = {

  fit_nodat <- orsf(formula = Surv(time, status) ~ .,
                    data = pbc_orsf,
                    n_tree = 1,
                    attach_data = FALSE)

  expect_error(
   orsf_pd_oob(fit_nodat, pred_spec = list(bili = c(0.8))),
   regexp = 'no data'
  )
 }
)

fit <- fit_standard_pbc$fast

test_that(
 "user cant supply empty pred_spec",
 code = {
  expect_error(
   orsf_ice_oob(fit, pred_spec = list()),
   regexp = 'pred_spec is empty'
  )
 }
)

test_that(
 "user cant supply pred_spec with non-matching names",
 code = {
  expect_error(
   orsf_ice_oob(fit,
                pred_spec = list(bili = 1:5,
                                 nope = c(1,2),
                                 no_sir = 1),
                pred_horizon = 1000),
   regexp = 'nope and no_sir'
  )
 }
)

bad_value_lower <- quantile(pbc_orsf$bili, probs = 0.01)
bad_value_upper <- quantile(pbc_orsf$bili, probs = 0.99)

test_that(
 "user cant supply pred_spec with values out of bounds",
 code = {
  expect_error(
   orsf_pd_new(fit,
               new_data = pbc_orsf,
               pred_spec = list(bili = c(bad_value_lower, 1:10, bad_value_upper)),
               pred_horizon = 1000),
   regexp = 'values for bili'
  )
 }
)

funs <- list(
 ice_new = orsf_ice_new,
 ice_inb = orsf_ice_inb,
 ice_oob = orsf_ice_oob,
 pd_new = orsf_pd_new,
 pd_inb = orsf_pd_inb,
 pd_oob = orsf_pd_oob
)

args_loop <- args_grid <- list(
 object = fit,
 pred_spec = list(bili = 1:4, sex = c("m", "f")),
 new_data = pbc_test,
 pred_horizon = 1000,
 pred_type = 'risk',
 na_action = 'fail',
 expand_grid = TRUE,
 prob_values = c(0.025, 0.50, 0.975),
 prob_labels = c("lwr", "medn", "upr"),
 boundary_checks = TRUE,
 n_thread = 3
)

args_loop$expand_grid <- FALSE

for(i in seq_along(funs)){

 f_name <- names(funs)[i]

 formals <- setdiff(names(formals(funs[[i]])), '...')

 for(pred_type in c('mort')){
 # for(pred_type in setdiff(pred_types_surv, c('leaf', 'mort'))){

  args_grid$pred_type = pred_type
  args_loop$pred_type = pred_type

  pd_object_grid <- do.call(funs[[i]], args = args_grid[formals])
  pd_object_loop <- do.call(funs[[i]], args = args_loop[formals])

  test_that(
   desc = paste('pred_spec data are returned on the original scale',
                ' for orsf_', f_name, sep = ''),
   code = {
    expect_equal(unique(pd_object_grid$bili), 1:4)
    expect_equal(unique(pd_object_loop[variable == 'bili', value]), 1:4)
   }
  )

  test_that(
   desc = paste(f_name, 'returns a data.table'),
   code = {
    expect_s3_class(pd_object_grid, 'data.table')
    expect_s3_class(pd_object_loop, 'data.table')
   }
  )

  test_that(
   desc = 'output is named correctly',
   code = {

    if(f_name %in% c("ice_new", "ice_inb", "ice_oob")){
     expect_true('id_variable' %in% names(pd_object_grid))
     expect_true('id_variable' %in% names(pd_object_loop))
     expect_true('id_row' %in% names(pd_object_grid))
     expect_true('id_row' %in% names(pd_object_loop))
    }

    expect_true('variable' %in% names(pd_object_loop))
    expect_true('value' %in% names(pd_object_loop))

    vars <- names(args_loop$pred_spec)
    expect_true(all(vars %in% names(pd_object_grid)))
    expect_true(all(vars %in% unique(pd_object_loop$variable)))

    if(pred_type == 'mort'){
     expect_false('pred_horizon' %in% names(pd_object_grid))
     expect_false('pred_horizon' %in% names(pd_object_loop))
    }

   }
  )



 }


}


pd_vals_ice <- orsf_ice_new(
 fit,
 new_data = pbc_orsf,
 pred_spec = list(bili = 1:4),
 pred_horizon = 1000
)

pd_vals_smry <- orsf_pd_new(
 fit,
 new_data = pbc_orsf,
 pred_spec = list(bili = 1:4),
 pred_horizon = 1000
)

test_that(
 'ice values summarized are the same as pd values',
 code = {

  pd_vals_check <- pd_vals_ice[, .(medn = median(pred)), by = id_variable]

  expect_equal(
   pd_vals_check$medn,
   pd_vals_smry$medn
  )

 }
)


test_that(
 'No missing values in output',
 code = {

  # expect_false(any(is.na(pd_vals_ice)))
  # expect_false(any(is.nan(as.matrix(pd_vals_ice))))
  # expect_false(any(is.infinite(as.matrix(pd_vals_ice))))

  expect_false(any(is.na(pd_vals_smry)))
  expect_false(any(is.nan(as.matrix(pd_vals_smry))))
  expect_false(any(is.infinite(as.matrix(pd_vals_smry))))
 }
)

test_that(
 'multi-valued horizon inputs are allowed',
 code = {

  pd_smry_multi_horiz <- orsf_pd_oob(
   fit,
   pred_spec = list(bili = 1),
   pred_horizon = c(1000, 2000, 3000)
  )

  # risk must increase or remain steady over time
  expect_lte(pd_smry_multi_horiz$mean[1], pd_smry_multi_horiz$mean[2])
  expect_lte(pd_smry_multi_horiz$mean[2], pd_smry_multi_horiz$mean[3])

  expect_lte(pd_smry_multi_horiz$medn[1], pd_smry_multi_horiz$medn[2])
  expect_lte(pd_smry_multi_horiz$medn[2], pd_smry_multi_horiz$medn[3])

  pd_ice_multi_horiz <- orsf_ice_oob(
   fit,
   pred_spec = list(bili = 1),
   pred_horizon = c(1000, 2000, 3000)
  )

  ice_check <- pd_ice_multi_horiz[, .(m = mean(pred, na.rm=TRUE)), by = pred_horizon]

  expect_equal(ice_check$m, pd_smry_multi_horiz$mean)

 }


)

# These tests are kept commented out and run locally
# I dont want to suggest pdp package in DESCRIPTION just for testing
# library(pdp)
#
# pred_aorsf <- function(object, newdata) {  # see ?predict.orsf_fit
#  as.numeric(predict(object, newdata, pred_horizon = 1000))
# }
#
# pd_reference <- partial(fit,
#                         pred.var = "bili",
#                         pred.grid = data.frame(bili = 1:5),
#                         pred.fun = pred_aorsf,
#                         plot = FALSE,
#                         ice = TRUE,
#                         train = pbc_orsf)
#
# pd_refsort <- as.data.table(pd_reference[order(pd_reference$bili), ])
#
# pd_refsort[, .(bili_mean = mean(yhat),
#                bili_median = median(yhat)),
#            by = bili]
#
# pbc_tmp <- pbc_orsf
#
# pbc_tmp$bili <- 1
#
# pred_bili_1 <- predict(fit, new_data = pbc_tmp, pred_horizon = 1000)
#
# pd_refsort[bili==1]
#
# head(pred_bili_1)
#
# pd_spec <- list(bili = 1:5)
#
# pd_bcj <- orsf_ice_inb(fit,
#                        pred_spec = pd_spec,
#                        pred_horizon = 1000,
#                        expand_grid = TRUE)
#
# pd_smry <- orsf_pd_inb(fit,
#                        pred_spec = pd_spec,
#                        pred_horizon = 1000)
#
# pd_bcj_check <- pd_bcj[, .(bili_mean = mean(pred),
#                            bili_median = median(pred)),
#                        by = bili]
#
# test_that(
#  desc = "pd_smry matches pd_ice",
#  code = {
#   expect_equal(pd_bcj_check$bili_mean, pd_smry$mean)
#   expect_equal(pd_bcj_check$bili_median, pd_smry$medn)
#  }
# )
#
# test_that(
#  "aorsf pd matches pdp package",
#  code = {
#   expect_equal(pd_bcj$pred, pd_refsort$yhat)
#  }
# )

# bili_seq <- seq(1, 5, length.out=20)
#
# microbenchmark::microbenchmark(
#  pd_reference = partial(fit,
#                         pred.var = "bili",
#                         pred.grid = data.frame(bili = bili_seq),
#                         pred.fun = pred_aorsf,
#                         plot = FALSE,
#                         ice = TRUE,
#                         train = pbc_orsf),
#  pd_aorsf = orsf_ice_new(fit,
#                          new_data = pbc_orsf,
#                          pred_spec = list(bili=bili_seq),
#                          pred_horizon = 1000,
#                          expand_grid = TRUE)
# )


