
object <- orsf(formula = Surv(time, status) ~ . - id,
               data_train = pbc_orsf,
               mtry = 5,
               n_split = 10,
               n_tree = 500,
               oobag_pred = T,
               oobag_time = 2500,
               leaf_min_obs = 10)

test_that(
 "user cant supply empty pd_spec",
 code = {
  expect_error(
   orsf_pd_ice(object,
               pd_data = pbc_orsf,
               pd_spec = list(),
               times = 1000,
               oobag = FALSE),
   regexp = 'pd_spec is empty'
  )
 }
)

test_that(
 "user cant supply pd_spec with non-matching names",
 code = {
  expect_error(
   orsf_pd_ice(object,
               pd_data = pbc_orsf,
               pd_spec = list(bili = 1:5,
                              nope = c(1,2),
                              no_sir = 1),
               times = 1000,
               oobag = FALSE),
   regexp = 'nope and no_sir'
  )
 }
)

bad_value_lower <- quantile(pbc_orsf$bili, probs = 0.01)
bad_value_upper <- quantile(pbc_orsf$bili, probs = 0.99)

test_that(
 "user cant supply pd_spec with values out of bounds",
 code = {
  expect_error(
   orsf_pd_ice(object,
               pd_data = pbc_orsf,
               pd_spec = list(bili = c(bad_value_lower, 1:10, bad_value_upper)),
               times = 1000,
               oobag = FALSE),
   regexp = 'values for bili'
  )
 }
)


# These tests are kept commented out and run locally
# I dont want to suggest pdp package in DESCRIPTION just for testing
# library(pdp)
#
# pred_aorsf <- function(object, newdata) {  # see ?predict.aorsf
#  as.numeric(predict(object, newdata, times = 1000))
# }
# pd_reference <- partial(object,
#                         pred.var = "bili",
#                         pred.grid = data.frame(bili = 1:5),
#                         pred.fun = pred_aorsf,
#                         plot = FALSE,
#                         ice = TRUE,
#                         train = pbc_orsf)
#
# pd_refsort <- pd_reference[order(pd_reference$bili), ]
#
# pd_spec <- list(bili = 1:5)
#
# pd_bcj <- orsf_pd_ice(object,
#                       pd_data = pbc_orsf,
#                       pd_spec = pd_spec,
#                       times = 1000,
#                       expand_grid = TRUE,
#                       oobag = FALSE)
#
# pd_smry <- orsf_pd_summary(object,
#                            pd_data = pbc_orsf,
#                            pd_spec = pd_spec,
#                            times = 1000,
#                            oobag = FALSE)
#
# pd_bcj_check <- pd_bcj[, .(bili_mean = mean(pred),
#                            bili_median = median(pred)),
#                        by = bili]
#
# test_that(
#  desc = "pd_smry matches pd_ice",
#  code = {
#   expect_equal(pd_bcj_check$bili_mean, pd_smry$mean)
#   expect_equal(pd_bcj_check$bili_median, pd_smry$median)
#  }
# )
#
# test_that(
#  "c fun matches R wrapper with pdp package",
#  code = {
#   expect_equal(pd_bcj$pred, pd_refsort$yhat)
#  }
# )



