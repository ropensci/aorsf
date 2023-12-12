
# survival ----

fit <- fit_standard_pbc$fast

pd_object_grid <- orsf_pd_oob(object = fit,
                              pred_spec = pred_spec_auto(sex, bili),
                              pred_horizon = c(1000, 2000))

pd_object_loop <- orsf_pd_oob(object = fit,
                              expand_grid = FALSE,
                              pred_spec = pred_spec_auto(sex, bili),
                              pred_horizon = c(1000, 2000))

test_that(
 desc = 'pred_spec data are returned on the original scale',
 code = {
  expect_equal(unique(pd_object_grid$bili),
               fit$get_var_bounds('bili'))

  expect_equal(unique(na.omit(pd_object_loop$value)),
               fit$get_var_bounds('bili'))
 }
)

test_that(
 desc = 'output is a data.table',
 code = {
  expect_s3_class(pd_object_grid, 'data.table')
  expect_s3_class(pd_object_loop, 'data.table')
 }
)


# classification ----

fit <- fit_standard_penguin_species$fast

pd_object_grid <- orsf_pd_new(object = fit,
                              new_data = penguins_orsf,
                              pred_spec = pred_spec_auto(bill_length_mm,
                                                         bill_depth_mm))

pd_object_loop <- orsf_pd_new(object = fit,
                              new_data = penguins_orsf,
                              expand_grid = FALSE,
                              pred_spec = pred_spec_auto(bill_length_mm,
                                                         bill_depth_mm))

test_that(
 desc = "probability values are bounded",
 code = {
  expect_true(all(pd_object_grid$mean <= 1))
  expect_true(all(pd_object_grid$mean >= 0))
  expect_true(all(pd_object_loop$mean <= 1))
  expect_true(all(pd_object_loop$mean >= 0))
 }
)


# regression ----

fit <- fit_standard_penguin_bills$fast

pd_object_grid <- orsf_pd_inb(object = fit,
                              pred_spec = pred_spec_auto(species, island),
                              pred_horizon = c(1000, 2000))

pd_object_loop <- orsf_pd_inb(object = fit,
                              expand_grid = FALSE,
                              pred_spec = pred_spec_auto(species, island),
                              pred_horizon = c(1000, 2000))

test_that(
 desc = "levels are converted to character values only if needed",
 code = {
  expect_false(is.character(pd_object_grid$species))
  expect_true(is.character(pd_object_loop$level))
 }
)


# general ----

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
 "cant supply empty pred_spec",
 code = {
  expect_error(
   orsf_ice_oob(fit, pred_spec = list()),
   regexp = 'pred_spec is empty'
  )
 }
)

test_that(
 "cant supply pred_spec with non-matching names",
 code = {
  expect_error(
   orsf_ice_oob(fit,
                pred_spec = list(bili = 1:5,
                                 nope = c(1,2),
                                 no_sir = 1),
                pred_horizon = 1000),
   regexp = 'nope and no_sir'
  )

  expect_error(
   orsf_ice_oob(fit,
                pred_spec = pred_spec_auto(bili, nope, no_sir),
                pred_horizon = 1000),
   regexp = 'nope and no_sir'
  )
 }
)

bad_value_lower <- quantile(pbc_orsf$bili, probs = 0.01)
bad_value_upper <- quantile(pbc_orsf$bili, probs = 0.99)

test_that(
 "cant supply pred_spec with values out of bounds",
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
 new_data = pbc_test,
 pred_spec = list(bili = 1:4),
 pred_horizon = 1000
)


pd_vals_smry <- orsf_pd_new(
 fit,
 new_data = pbc_test,
 pred_spec = list(bili = 1:4),
 pred_horizon = 1000
)

test_that(
 'ice values summarized are the same as pd values',
 code = {

  grps <- split(pd_vals_ice, pd_vals_ice$id_variable)
  pd_vals_check <- sapply(grps, function(x) median(x$pred))

  expect_equal(
   as.numeric(pd_vals_check),
   pd_vals_smry$medn
  )

 }
)



test_that(
 'No missing values in summary output',
 code = {
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
   pred_type = 'risk',
   pred_spec = list(bili = 1),
   pred_horizon = c(1000, 2000, 3000)
  )

  # risk monotonically increases
  expect_lte(pd_smry_multi_horiz$mean[1], pd_smry_multi_horiz$mean[2])
  expect_lte(pd_smry_multi_horiz$mean[2], pd_smry_multi_horiz$mean[3])
  expect_lte(pd_smry_multi_horiz$medn[1], pd_smry_multi_horiz$medn[2])
  expect_lte(pd_smry_multi_horiz$medn[2], pd_smry_multi_horiz$medn[3])

  pd_ice_multi_horiz <- orsf_ice_oob(
   fit,
   pred_type = 'risk',
   pred_spec = list(bili = 1),
   pred_horizon = c(1000, 2000, 3000)
  )

  grps <- split(pd_ice_multi_horiz, pd_ice_multi_horiz$pred_horizon)

  ice_check <- sapply(grps, function(x) mean(x$pred, na.rm=TRUE))

  expect_equal(as.numeric(ice_check), pd_smry_multi_horiz$mean)

 }


)

# # These tests are kept commented out and run locally
# # I dont want to suggest pdp package in DESCRIPTION just for testing
# library(pdp)
#
# pred_aorsf <- function(object, newdata) {  # see ?predict.ObliqueForest
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
#
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


