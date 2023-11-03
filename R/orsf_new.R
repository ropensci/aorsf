
orsf_new <- function(...){

 # input checks ----

 inputs <- list(...)

 # this makes code below a little less tedious
 data <- inputs$data

 do.call(check_orsf_inputs, args = inputs)

 # untrained at initialization
 inputs$trained <- FALSE

 # Y attributes ----

 names_y_data <- check_var_names(data, all.vars(inputs$formula[[2]]))

 if(any(is.na(select_cols(data, names_y_data))))
  stop("Please remove missing values from the outcome variable(s)",
       call. = FALSE)

 vt = c('numeric', 'integer', 'units', 'factor')

 types_y_data <- check_var_types(data, names_y_data, valid_types = vt)

 unit_y_names <- names_y_data[types_y_data == 'units']

 ui_y <- unit_info(data = data, .names = unit_y_names)

 # X attributes ----

 # info about formula variables
 formula_terms <- suppressWarnings(
  stats::terms(x = inputs$formula, data = data)
 )

 if(attr(formula_terms, 'response') == 0)
  stop("formula should have a response", call. = FALSE)

 # store x names and types from training data
 names_x_data <- check_var_names(data, attr(formula_terms, 'term.labels'))

 inputs$rows_x_cc <- stats::complete.cases(select_cols(data, names_x_data))

 vt <- c('numeric', 'integer', 'units', 'factor', 'ordered')

 types_x_data <- check_var_types(data, names_x_data, valid_types = vt)

 # store the factor information
 fi <- fctr_info(data, names_x_data)

 # store the number of columns in X matrix
 # X matrix isn't created yet but we can infer the number of cols in X
 # by using factor info and the length of numeric columns.

 n_col_x_nominal <- sum(
  vapply(fi$lvls, function(x) length(x) - 1L, integer(1))
 )

 n_col_x_numeric <- sum(
  grepl(pattern = "^integer$|^numeric$|^units$", x = types_x_data)
 )

 inputs$n_col_x <- n_col_x_nominal + n_col_x_numeric

 unit_x_names <- names_x_data[types_x_data == 'units']

 ui_x <- unit_info(data = data, .names = unit_x_names)

 # X & Y attributes ----
 names_all <- c(names_x_data, names_y_data)

 # check for infinite or all missing predictors
 check_var_values(data, names_all)

 # check factors in x data
 fctr_check(data, names_all)
 fctr_id_check(data, names_all)

 # store attributes ----
 inputs$types_x_data <- types_x_data
 inputs$names_x_data <- names_x_data
 inputs$types_y_data <- types_y_data
 inputs$names_y_data <- names_y_data
 inputs$fctr_info <- fi
 inputs$unit_info <- c(ui_y, ui_x)

 if(is.null(inputs$mtry)){
  inputs$mtry <-  ceiling(sqrt(inputs$n_col_x))
 }

 inputs$n_obs <- ifelse(test = inputs$na_action == 'omit',
                        yes  = length(inputs$rows_x_cc),
                        no   = nrow(data))

 # delayed checks ----

 check_arg_lteq(
  arg_value = inputs$mtry,
  arg_name = 'mtry',
  bound = inputs$n_col_x,
  append_to_msg = "(number of columns in the one-hot encoded x-matrix)"
 )

 if(is.null(inputs$control$lincomb_df_target)){

  inputs$control$lincomb_df_target <- inputs$mtry

 } else {

  check_arg_lteq(
   arg_value = inputs$control$lincomb_df_target,
   arg_name = 'df_target',
   bound = inputs$mtry,
   append_to_msg = "(number of randomly selected predictors)"
  )

 }


 # set weights as 1 if user did not supply them.
 # length of weights depends on how missing are handled.
 if(is.null(inputs$weights)){

  inputs$weights <- rep(1, inputs$n_obs)

 } else {

  check_arg_length(
   arg_value = inputs$weights,
   arg_name  = 'weights',
   expected_length = inputs$n_obs
  )

 }

 check_arg_lteq(arg_value = inputs$leaf_min_obs,
                arg_name = "leaf_min_obs",
                bound = round(inputs$n_obs / 2),
                append_to_msg = "(number of observations divided by 2)")

 check_arg_lt(arg_value = inputs$split_min_obs,
              arg_name = "split_min_obs",
              bound = inputs$n_obs,
              append_to_msg = "(number of observations)")

 # tree attributes ----

 inputs$tree_type = infer_tree_type(names_y_data, data)

 if(is.null(inputs$tree_seeds)){

  inputs$tree_seeds <- sample(1e6, size = 1)

 }

 if(length(inputs$tree_seeds) == 1){

  if(inputs$n_tree > 1) set.seed(inputs$tree_seeds)

  inputs$tree_seeds <- sample(inputs$n_tree*10, size = inputs$n_tree)

 }

 # out-of-bag attributes ----

 inputs$oobag_pred <- inputs$oobag_pred_type != "none"

 if(is.null(inputs$oobag_fun)){

  inputs$f_oobag_eval <- function(x) x
  inputs$type_oobag_eval <- 'none'

  if(inputs$oobag_pred){

   inputs$type_oobag_eval <- switch(inputs$tree_type,
                                    'survival' = 'cstat',
                                    'classification' = 'cstat',
                                    'regression' = 'rsq')
  }

 } else {

  inputs$type_oobag_eval <- 'user'
  names(inputs)[names(inputs) == 'oobag_fun'] <- 'f_oobag_eval'

 }

 # can't evaluate the oobag predictions if they aren't aggregated
 if(inputs$oobag_pred_type == 'leaf'){
  inputs$type_oobag_eval <- 'none'
 }

 # transfer data to its own list after checking
 fields <- list(data = inputs$data)
 fields$forest = list()
 fields$importance = double(length = 0)
 fields$pred_oobag = matrix(ncol = 0, nrow = 0)
 fields$eval_oobag = list(
  stat_values = matrix(ncol = 0, nrow = 0),
  stat_type = switch(inputs$type_oobag_eval,
                     'none' = "None",
                     'cstat' = "C statistic",
                     'user' = "User-specified function")
 )

 # initialize new orsf_fit ----

 out <- structure(.Data = fields)

 # replace data with names in attribute list
 inputs[[1]] <- names(out)
 names(inputs)[1] <- "names"
 # attach all inputs as attributes
 attributes(out) <- inputs
 # indicate S3 class
 class(out) <- c(paste("orsf_fit", inputs$tree_type, sep = "_"),
                 "orsf_fit")

 # add tree-specific attributes ----

 internal_attributes <- orsf_new_internals(out)

 for(i in seq_along(internal_attributes)){
  attr(out, names(internal_attributes)[i]) <- internal_attributes[[i]]
 }

 out

}

orsf_new_internals <- function(object){
 UseMethod("orsf_new_internals", object)
}

orsf_new_internals.orsf_fit_survival <- function(object){

 y <- select_cols(object$data, get_names_y(object))

 if(get_na_action(object) == 'na_omit'){
  y <- y[get_rows_x_cc(object), ]
 }

 n_events <- collapse::fsum(y[, 2])

 check_arg_lteq(
  arg_value = get_leaf_min_events(object),
  arg_name = 'leaf_min_events',
  bound = round(n_events / 2),
  append_to_msg = "(number of events divided by 2)"
 )

 check_arg_lt(
  arg_value = get_split_min_events(object),
  arg_name = "split_min_events",
  bound = n_events,
  append_to_msg = "(number of events)"
 )

 # if pred_horizon is unspecified, provide sensible default
 ph <- get_oobag_pred_horizon(object) %||% collapse::fmedian(y[, 1])

 # order observations by event time and event status
 sorted <- collapse::radixorder(y[, 1], -y[, 2])

 max_time <- y[last_value(sorted), 1]

 list(
  oobag_pred_horizon = ph,
  sorted = sorted,
  max_time = max_time
 )

}
orsf_new_internals.orsf_fit_classification <- function(object){
 NULL
}
orsf_new_internals.orsf_fit_regression <- function(object){
 NULL
}

orsf_data_summarize <- function(object){

 names_x_data <- get_names_x(object)
 names_y_data <- get_names_y(object)
 weights      <- get_weights(object)
 data         <- object$data

 names_x_numeric <- grep(pattern = "^integer$|^numeric$|^units$",
                         x = get_types_x(object))

 means <- standard_deviations <- modes <- numeric_bounds <- NULL

 numeric_cols <- names_x_data[names_x_numeric]
 nominal_cols <- get_fctr_info(object)$cols

 if(!is_empty(nominal_cols)){

  modes <- vapply(
   select_cols(object$data, nominal_cols),
   collapse::fmode,
   FUN.VALUE = integer(1),
   w = weights
  )

 }

 if(!is_empty(numeric_cols)){

  numeric_data <- select_cols(data, numeric_cols)

  numeric_bounds <- matrix(
   data = c(
    collapse::fnth(numeric_data, 0.1, w = weights),
    collapse::fnth(numeric_data, 0.25, w = weights),
    collapse::fnth(numeric_data, 0.5, w = weights),
    collapse::fnth(numeric_data, 0.75, w = weights),
    collapse::fnth(numeric_data, 0.9, w = weights)
   ),
   nrow = 5,
   byrow = TRUE,
   dimnames = list(c('10%', '25%', '50%', '75%', '90%'),
                   names(numeric_data))
  )

  means <- collapse::fmean(numeric_data, w = weights)

  standard_deviations <- collapse::fsd(numeric_data, w = weights)

 }

 list(means = means,
      modes = modes,
      standard_deviations = standard_deviations,
      numeric_bounds = numeric_bounds)

}

orsf_prep_xyw <- function(object){

 data <- object$data
 names_x_data <- get_names_x(object)
 names_y_data <- get_names_y(object)

 if(any(is.na(select_cols(data, names_x_data)))){

  switch(
   get_na_action(object),

   'fail' = {
    stop("Please remove missing values from data, or impute them.",
         call. = FALSE)
   },

   'omit' = {
    data <- data[get_rows_x_cc(object), ]
   },

   'impute_meanmode' = {

    data <- data_impute(data,
                        cols = names_x_data,
                        values = c(as.list(get_means(object)),
                                   as.list(get_modes(object))))
   }

  )

 }

 # dont use generic here because data may have been modified
 y <- switch(
  get_tree_type(object),
  'survival' = prep_y_surv(data, names_y_data),
  'classification' = prep_y_clsf(data, names_y_data)
 )

 x <- prep_x(data,
             get_fctr_info(object),
             names_x_data,
             get_means(object),
             get_standard_deviations(object))

 list(x = x, y = y, w = get_weights(object))

}


orsf_clean_fit <- function(object){

 importance_clean <- NULL
 pred_oobag_clean <- NULL
 stat_values_clean <- NULL

 if(get_importance(object) != "none"){

  importance_clean <- object$importance

  rownames(importance_clean) <- get_names_x(object, ref_code_names = TRUE)

  importance_clean <-
   rev(importance_clean[order(importance_clean), , drop=TRUE])

 }

 if(get_oobag_pred(object)){

  pred_oobag_clean <- object$pred_oobag

  if(get_oobag_pred_type(object) == 'leaf'){

   all_rows <- seq(get_n_obs(object))

   for(i in seq(get_n_tree(object))){

    rows_inbag <- setdiff(all_rows, object$forest$rows_oobag[[i]]+1)
    pred_oobag_clean[rows_inbag, i] <- NA

   }

  }

  pred_oobag_clean[is.nan(pred_oobag_clean)] <- NA_real_

  if(get_tree_type(object) == 'survival'){

   # put the oob predictions into the same order as the training data.
   unsorted <- collapse::radixorder(get_sorted(object))
   pred_oobag_clean <- pred_oobag_clean[unsorted, , drop = FALSE]

   # mortality predictions should always be 1 column
   # b/c they do not depend on the prediction horizon
   if(get_oobag_pred_type(object) == 'mort'){

    stat_values_clean <- object$eval_oobag$stat_values
    stat_values_clean <- stat_values_clean[, 1L, drop = FALSE]
    pred_oobag_clean <- pred_oobag_clean[, 1L, drop = FALSE]

   }

  }


 }

 out <- list(importance = importance_clean,
             pred_oobag = pred_oobag_clean,
             stat_values = stat_values_clean)

 internals <- orsf_clean_fit_internal(object)

 c(out, internals)


}

orsf_clean_fit_internal <- function(object){
 UseMethod("orsf_clean_fit_internal", object)
}

orsf_clean_fit_internal.orsf_fit_survival <- function(object){

 list(pred_horizon = get_oobag_pred_horizon(object))

}

orsf_clean_fit_internal.orsf_fit_classification <- function(object){
 NULL
}

orsf_clean_fit_internal.orsf_fit_regression <- function(object){
 NULL
}


