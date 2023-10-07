#' @srrstats {G5.2} *Appropriate warnings/error explicitly demonstrated through tests.*
#' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*
#' @srrstats {G5.3} *Test that fits returned contain no missing (`NA`) or undefined (`NaN`, `Inf`) values.*
#' @srrstats {G5.8, G5.8d} **Edge condition tests** * an error is thrown when partial dependence functions are asked to predict estimates outside of boundaries determined by the aorsf model's training data*

test_that(
 desc = "oob stops if there are no data",
 code = {

  fit_nodat <- orsf(formula = Surv(time, status) ~ .,
                    data = pbc_orsf,
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
 'pred_spec data are returned on the original scale',

 code = {

  expect_equal(unique(pd_vals_ice$bili), 1:4)
  expect_equal(unique(pd_vals_smry$bili), 1:4)

 }

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

  expect_false(any(is.na(pd_vals_ice)))
  expect_false(any(is.nan(as.matrix(pd_vals_ice))))
  expect_false(any(is.infinite(as.matrix(pd_vals_ice))))

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

  ice_check <- pd_ice_multi_horiz[, .(m = mean(pred)), by = pred_horizon]

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



