
orsf_new <- function(...){

 inputs <- list(...)

 # this makes code below a little less tedious
 data <- inputs$data

 # number of observations total and complete case
 inputs$n_obs <- nrow(data)
 inputs$n_obs_cc <- nrow(collapse::na_omit(data))

 do.call(check_orsf_inputs, args = inputs)

 # untrained at initialization
 inputs$trained <- FALSE

 # info about formula variables
 formula_terms <- suppressWarnings(
  stats::terms(x = inputs$formula, data = data)
 )

 # supervised forests only
 if(attr(formula_terms, 'response') == 0)
  stop("formula should have a response", call. = FALSE)

 # store x names and types from training data
 names_x_data <- check_var_names(data, attr(formula_terms, 'term.labels'))

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

 if(is.null(inputs$mtry)){
  inputs$mtry <-  ceiling(sqrt(ncol(inputs$n_col_x)))
 }

 # run additional checks that are dependent on the data
 check_arg_lteq(
  arg_value = inputs$mtry,
  arg_name = 'mtry',
  bound = inputs$n_col_x,
  append_to_msg = "(number of columns in the one-hot encoded x-matrix)"
 )

 if(is.null(inputs$control$lincomb_df_target)){

  inputs$control$lincomb_df_target <- mtry

 } else {

  check_arg_lteq(
   arg_value = inputs$control$lincomb_df_target,
   arg_name = 'df_target',
   bound = mtry,
   append_to_msg = "(number of randomly selected predictors)"
  )

 }


 unit_x_names <- names_x_data[types_x_data == 'units']

 ui_x <- unit_info(data = data, .names = unit_x_names)


 # store y names and types from training data
 names_y_data <- check_var_names(data, all.vars(inputs$formula[[2]]))

 if(any(is.na(select_cols(data, names_y_data))))
  stop("Please remove missing values from the outcome variable(s)",
       call. = FALSE)

 vt = c('numeric', 'integer', 'units', 'factor')

 types_y_data <- check_var_types(data, names_y_data, valid_types = vt)

 unit_y_names <- names_y_data[types_y_data == 'units']

 ui_y <- unit_info(data = data, .names = unit_y_names)

 # checks for both x and y
 names_all <- c(names_x_data, names_y_data)

 # check for infinite or all missing predictors
 check_var_values(data, names_all)

 # check factors in x data
 fctr_check(data, names_all)
 fctr_id_check(data, names_all)

 # store x attributes
 inputs$types_x_data <- types_x_data
 inputs$names_x_data <- names_x_data
 # store y attributes
 inputs$types_y_data <- types_y_data
 inputs$names_y_data <- names_y_data
 # store shared attributes
 inputs$fctr_info <- fi
 inputs$unit_info <- c(ui_y, ui_x)

 # determine and store the tree type
 inputs$tree_type = infer_tree_type(names_y_data, data)

 # out-of-bag attributes
 inputs$oobag_pred <- inputs$oobag_pred_type != "none"

 if(is.null(inputs$oobag_fun)){
  inputs$type_oobag_eval <- if(inputs$oobag_pred) 'cstat' else 'none'
  inputs$f_oobag_eval <- function(x) x
 } else {
  inputs$type_oobag_eval <- 'user'
  names(inputs)[names(inputs) == 'oobag_fun'] <- 'f_oobag_eval'
 }

 # can't evaluate the oobag predictions if they aren't aggregated
 if(inputs$oobag_pred_type == 'leaf'){
  inputs$type_oobag_eval <- 'none'
 }

 # rename weights so as not to confuse with other weights
 names(inputs)[names(inputs) == 'oobag_fun'] <- 'f_oobag_eval'

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

 out <- structure(.Data = fields)

 # replace data with names in attribute list
 inputs[[1]] <- names(out)
 names(inputs)[1] <- "names"
 # attach all inputs as attributes
 attributes(out) <- inputs

 class(out) <- c(paste("orsf_fit", inputs$tree_type, sep = "_"),
                 "orsf_fit")

 out

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

orsf_prep <- function(){

}
