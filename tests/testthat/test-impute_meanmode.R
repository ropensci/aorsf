
suppressPackageStartupMessages({
 library(collapse)
})

pbc_miss <- as.data.table(survival::pbc) %>%
 .[status > 0, status := status - 1] %>%
 .[, stage := factor(stage, ordered = TRUE)] %>%
 .[, trt := factor(trt,
                   levels = c(1, 2),
                   labels = c('d_penicill_main',
                              'placebo'))] %>%
 ftransformv(vars = c(ascites, hepato, spiders, edema), FUN = factor)

fit_miss <- orsf(pbc_miss,
                 time + status ~ . - id,
                 na_action = 'impute_meanmode')

expect_equal(
 sum(complete.cases(fit_miss$data)),
 nrow(pbc_miss)
)

impute_values <- c(as.list(get_means(fit_miss)),
                   as.list(get_modes(fit_miss)))

test_that(
 desc = "imputation does not coerce columns to new types",
 code = {
  for(i in names(pbc_miss)){
   expect_equal(
    typeof(pbc_miss[[i]]),
    typeof(fit_miss$data[[i]])
   )
  }
 }
)

test_that(
 desc = "integer cols imputed by coercing imputed value to integer",
 code = {
  chol_na <- whichNA(pbc_miss$chol)
  expect_equal(
   fit_miss$data$chol[chol_na],
   rep(as.integer(impute_values$chol), length(chol_na))
  )
 }
)

test_that(
 desc = "factor cols imputed with the level corresponding to stored int",
 code = {
  trt_na <- whichNA(pbc_miss$trt)
  expect_true(
   all(fit_miss$data$trt[trt_na] == levels(pbc_miss$trt)[impute_values$trt])
  )
 }
)

test_that(
 desc = "doubles imputed with the stored double",
 code = {
  alk_na <- whichNA(pbc_miss$alk.phos)
  expect_true(
   all(fit_miss$data$alk.phos[alk_na] == impute_values$alk.phos)
  )
 }
)

test_that(
 desc = 'dimensions of predicted output match expectations',
 code = {

  pred_omit <- predict(fit_miss,
                       new_data = pbc_miss,
                       na_action = 'omit')

  expect_equal(
   nrow(pred_omit),
   sum(complete.cases(pbc_miss))
  )

  pred_pass <- predict(fit_miss,
                       new_data = pbc_miss,
                       na_action = 'pass')

  expect_equal(
   complete.cases(pred_pass),
   complete.cases(pbc_miss)
  )

  pred_impute <- predict(fit_miss,
                         new_data = pbc_miss,
                         na_action = 'impute_meanmode')

  expect_equal(
   sum(complete.cases(pred_impute)),
   nrow(pbc_miss)
  )

 }
)




#
# data <- survival::pbc
#
# numeric_cols <- names(which(vapply(data, is.numeric, logical(1))))
# categorical_cols <- names(which(vapply(data, is.factor, logical(1))))
#
# impute_values <- c(get_means(data, numeric_cols),
#                    get_modes(data, categorical_cols))
#
# for(i in names(impute_values)){
#
#  if(i %in% numeric_cols)
#   expect_equal(impute_values[[i]], mean(data[[i]], na.rm = TRUE))
#
#  if(i %in% categorical_cols)
#   expect_equal(impute_values[[i]], mode(data[[i]]))
#
# }
#
# missing_index <- is.na(data)
# colnames(missing_index) <- colnames(data)
#
# data_imputed <- impute_meanmode(data,
#                                 cols = names(data),
#                                 values = impute_values)
#
# impute_values_collapse <- impute_values
# impute_values_collapse$sex <- 2L
#
# test_that(
#  desc = "means and modes are placed into data after imputation",
#  code = {
#   for(col in names(data)){
#
#    if( any(missing_index[, col]) ){
#
#     rows_imputed <- which(missing_index[, col])
#
#     expect_equal(
#      data_imputed[rows_imputed, col],
#      rep(impute_values[[col]], length(rows_imputed))
#     )
#    }
#   }
#  }
# )
#
