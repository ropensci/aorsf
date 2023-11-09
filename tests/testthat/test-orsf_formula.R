

#' @srrstats {G5.2} *Appropriate error behaviour is explicitly demonstrated through tests.*
#' @srrstats {G5.2b} *Tests demonstrate conditions which trigger error messages.*
test_that(
 desc = 'formula inputs are vetted',
 code = {

  # set id to a factor so that it can trigger the id error
  pbc_orsf$id <- factor(pbc_orsf$id)

  expect_error(orsf(pbc_orsf, Surv(time, status) ~ unknown_variable + bili),
               'not found in data')

  expect_error(orsf(pbc_orsf, Surv(time, status) ~ bili + factor(hepato)),
               'unrecognized')

  expect_error(orsf(pbc_orsf, Surv(time, status) ~ bili * ascites),
               'unrecognized')

  expect_error(orsf(pbc_orsf, Surv(time, status) ~ bili + id),
               'id variable?')

  expect_error(orsf(pbc_orsf, Surv(time, not_right) ~ . - id),
               'not_right')

  expect_error(orsf(pbc_orsf, Surv(not_right, status) ~ . - id),
               'not_right')

  expect_error(orsf(pbc_orsf, Surv(start, time, status) ~ .),
               'should have at most two variables')

  expect_error(orsf(pbc_orsf, Surv(time, id) ~ . -id),
               'detected >1 event type')

  expect_error(orsf(pbc_orsf, ~ .), 'must be two sided')

  expect_error(orsf(pbc_orsf, time + status ~ id + bili),
               "as many levels as there are rows")

 }
)

test_that(
 desc = 'long formulas with repetition are allowed',
 code = {

  x_vars <- c(setdiff(names(pbc_orsf), c('time', 'status', 'id')))

  long_rhs <- paste(x_vars, collapse = ' + ')

  long_rhs <- rep(long_rhs, 15)

  long_rhs <- paste(long_rhs, collapse = ' + ')

  f_long <- as.formula(paste("time + status ~", long_rhs))

  for(i in seq_along(fit_standard_pbc)){

   fit_long <- orsf(pbc_orsf,
                    formula = f_long,
                    control = controls[[i]],
                    n_tree = n_tree_test,
                    tree_seeds = seeds_standard)

   # fits the orsf as expected
   expect_s3_class(fit_long, 'ObliqueForest')
   # keeps unique names
   expect_equal(x_vars, fit_long$get_names_x())
   # is the same forest as standard
   expect_equal_leaf_summary(fit_long, fit_standard_pbc[[i]])

  }

 }
)

test_that(
 desc = "Surv objects in formula are used correctly",
 code = {

  pbc_surv <- Surv(pbc_orsf$time, pbc_orsf$status)

  pbc_surv_data <- cbind(pbc_orsf, surv_object = pbc_surv)

  for(i in seq_along(controls)){

   fit_surv <- orsf(
    pbc_surv_data,
    formula = surv_object ~ . - id - time - status,
    n_tree = n_tree_test,
    control = controls[[i]],
    tree_seed = seeds_standard
   )

   # name of surv object is correctly stored, values can be reproduced
   expect_equal(
    pbc_surv_data[[fit_surv$get_names_y()]],
    pbc_surv
   )

   # different formula but same as standard forest
   expect_equal_leaf_summary(fit_surv, fit_standard_pbc[[i]])
  }

 }
)

# test_that(
#  desc = "Status can be 0/1 or 1/2, or generally x/x+1",
#  code = {
#   for(i in seq(1:5)){
#
#   pbc_orsf$status <- pbc_orsf$status+1
#
#   for(j in seq_along(fit_standard_pbc)){
#
#    fit_status_modified <- orsf(pbc_orsf,
#                                time + status ~ . - id,
#                                n_tree = n_tree_test,
#                                control = controls[[j]],
#                                tree_seeds = seeds_standard)
#
#    expect_equal_leaf_summary(fit_status_modified, fit_standard_pbc[[j]])
#
#   }
#
#   expect_error(
#    orsf(pbc_orsf, Surv(status, time) ~ . - id),
#    'Did you enter'
#   )
#
#   }
#  }
# )
