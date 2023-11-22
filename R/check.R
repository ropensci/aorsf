
#' strict checks for inputs
#'
#' @param arg_value the object that is to be checked
#' @param arg_name the name of the object (used for possible error message)
#' @param expected_type what type of object should this be?
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd

check_arg_type <- function(arg_value, arg_name, expected_type){

 if('numeric' %in% expected_type)
  expected_type <- c(setdiff(expected_type, 'numeric'),
                     'double', 'integer')

 arg_type <- typeof(arg_value)

 type_match <-
  arg_type %in% expected_type | inherits(arg_value, expected_type)

 if (!type_match) {

  expected_types <- paste_collapse(x = expected_type,
                                   sep = ', ',
                                   last = ' or ')

  error_msg <-
   paste0(arg_name, " should have type <", expected_types, ">",
          " but instead has type <", arg_type, ">")

  stop(error_msg, call. = FALSE)

 }

}

#' strict checks for inputs
#'
#' @param arg_value the object that is to be checked
#' @param arg_name the name of the object (used for possible error message)
#' @param expected_uni what unique values should `arg_value` have?
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd

check_arg_uni <- function(uni, arg_name, expected_uni){

 # expected_in_uni <- all(expected_uni %in% uni)
 uni_in_expected <- all(uni %in% expected_uni)

 expected_values <- paste_collapse(x = expected_uni,
                                   sep = ', ',
                                   last = ' and ')

 if(!uni_in_expected){

  invalid_values <- paste_collapse(x = setdiff(uni, expected_uni),
                                   sep = ', ',
                                   last = ' and ')

  error_msg <-
   paste0(arg_name, " should contain values of ", expected_values,
          " but has values of ", invalid_values)

  stop(error_msg, call. = FALSE)


 }

}


#' strict checks for inputs
#'
#' @param arg_value the object that is to be checked
#' @param arg_name the name of the object (used for possible error message)
#' @param expected_type what length should `arg_value` have?
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
check_arg_length <- function(arg_value, arg_name, expected_length){

 arg_length <- length(arg_value)

 length_match <- arg_length %in% expected_length

 if (!length_match) {

  expected_lengths <- paste_collapse(x = expected_length,
                                     sep = ', ',
                                     last = ' or ')

  error_msg <-
   paste0(arg_name, " should have length <", expected_lengths, ">",
          " but instead has length <", arg_length, ">")

  stop(error_msg, call. = FALSE)

 }

}

#' strict checks for inputs
#'
#' @param arg_value the object that is to be checked
#' @param arg_name the name of the object (used for possible error message)
#' @param bound what bounds to use for `arg_value`?
#' @param relational_operator <, <=, >, or >=. The operator determines
#'   how bounds are checked.
#' @param append_to_msg a note to be added to the error message.
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
check_arg_bound <- function(arg_value, arg_name, bound,
                            relational_operator,
                            append_to_msg){

 .op <- switch(relational_operator,
               'gt' = `>`,
               'lt' = `<`,
               'gteq' = `>=`,
               'lteq' = `<=`)

 .lab <- switch(relational_operator,
                'gt' = ">",
                'lt' = "<",
                'gteq' = ">=",
                'lteq' = "<=")

 .neg <- switch(relational_operator,
                'gt' = "<=",
                'lt' = ">=",
                'gteq' = "<",
                'lteq' = ">")

 fails <- !.op(arg_value, bound)

 if(any(fails)){

  if(length(arg_value) == 1){

   error_msg <-
    paste0(arg_name, " = ", arg_value, " should be ", .lab, " ", bound)

  } else {

   first_offense <- min(which(fails))

   error_msg <- paste0(arg_name, " should be ", .lab, " ", bound, " but has",
                       " at least one value that is ", .neg, " ", bound,
                       " (see ", arg_name, "[", first_offense, "])")
  }

  if(!is.null(append_to_msg)){

   error_msg <- paste(error_msg, append_to_msg)

  }

  stop(error_msg, call. = FALSE)

 }

}

#' strict checks for inputs
#'
#'  argument is strictly greater than a bound.
#'
#' @inheritParams check_arg_bound
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
check_arg_gt <- function(arg_value, arg_name, bound, append_to_msg = NULL){
 check_arg_bound(arg_value,
                 arg_name,
                 bound,
                 relational_operator = 'gt',
                 append_to_msg)
}

#' strict checks for inputs
#'
#'  argument is strictly less than a bound.
#'
#' @inheritParams check_arg_bound
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
check_arg_lt <- function(arg_value, arg_name, bound, append_to_msg = NULL){
 check_arg_bound(arg_value,
                 arg_name,
                 bound,
                 relational_operator = 'lt',
                 append_to_msg)
}

#' strict checks for inputs
#'
#'  argument is greater than or equal to a bound.
#'
#' @inheritParams check_arg_bound
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
check_arg_gteq <- function(arg_value, arg_name, bound, append_to_msg = NULL){
 check_arg_bound(arg_value,
                 arg_name,
                 bound,
                 relational_operator = 'gteq',
                 append_to_msg)
}

#' strict checks for inputs
#'
#'  argument is less than or equal to a bound.
#'
#' @inheritParams check_arg_bound
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd

check_arg_lteq <- function(arg_value, arg_name, bound, append_to_msg = NULL){
 check_arg_bound(arg_value,
                 arg_name,
                 bound,
                 relational_operator = 'lteq',
                 append_to_msg)
}

#' strict checks for inputs
#'
#' @param arg_value the object that is to be checked
#' @param arg_name the name of the object (used for possible error message)
#' @param valid_options what are the valid inputs for `arg_value`?
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
check_arg_is_valid <- function(arg_value, arg_name, valid_options) {

 valid_arg <- arg_value %in% valid_options

 if (!valid_arg) {

  expected_values <- paste_collapse(x = valid_options,
                                    sep = ', ',
                                    last = ' or ')

  arg_values <- paste_collapse(x = arg_value,
                               sep = ', ',
                               last = ' or ')

  error_msg <- paste0(
   arg_name, " should be <", expected_values, ">",
   " but is instead <", arg_values, ">"
  )

  stop(error_msg, call. = FALSE)

 }

}

#' strict checks for inputs
#'
#' make sure the user has supplied an integer valued input.
#'
#' @param arg_value the object that is to be checked
#' @param arg_name the name of the object (used for possible error message)
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
check_arg_is_integer <- function(arg_name, arg_value){

 is_integer <- all(as.integer(arg_value) == arg_value)

 if(!is_integer){

  if(length(arg_value) == 1){
   error_msg <- paste0(arg_name, " should be an integer value",
                       " but instead has a value of ", arg_value)
  } else {

   first_offense <- min(which(as.integer(arg_value) != arg_value))

   error_msg <- paste0(arg_name, " should contain only integer values",
                       " but has at least one double value",
                       " (see ", arg_name, "[", first_offense, "])")

  }

  stop(error_msg, call. = FALSE)

 }

}

#' strict checks for inputs
#'
#' @param arg_value the object that is to be checked
#' @param arg_name the name of the object (used for possible error message)
#' @param expected_class what class should `arg_value` inherit from?
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
check_arg_is <- function(arg_value, arg_name, expected_class){

 arg_is <- inherits(arg_value, expected_class)

 if (!arg_is) {

  expected_classes <- paste_collapse(x = expected_class,
                                     sep = ', ',
                                     last = ' or ')

  arg_classes <- paste_collapse(x = class(arg_value),
                                sep = ', ',
                                last = ' or ')

  error_msg <- paste0(
   arg_name, " should inherit from class <", expected_classes, ">",
   " but instead inherits from <", arg_classes, ">"
  )

  stop(error_msg, call. = FALSE)

 }

}


#' check mis-typed arguments
#'
#' @param .dots ... from a call to .f
#' @param .f the function being called
#'
#' @return an error if you have mis-typed an arg
#' @noRd

check_dots <- function(.dots, .f){

 if(!is_empty(.dots)){

  .args <- setdiff(names(formals(.f)), '...')
  .dots <- names(.dots)

  for(i in seq_along(.dots)){

   .match_indices <- utils::adist(x = .dots[i],
                                  y = .args,
                                  fixed = TRUE,
                                  costs = c(ins = 1,
                                            del = 1,
                                            sub = 2))

   .match_index <- which.min(.match_indices)

   .dots[i] <- paste('  ', .dots[i],
                     ' is unrecognized - did you mean ',
                     .args[.match_index], '?', sep = '')
  }


  stop("there were unrecognized arguments:\n",
       paste(.dots, collapse = '\n'),
       call. = FALSE)

 }

}



#' Check variable types
#'
#' orsf() should only be run with certain types of variables. This function
#'   checks input data to make sure all variables have a primary (i.e., first)
#'   class that is within the list of valid options.
#'
#' @param data data frame with variables to be checked
#' @param .names names of variables in `data` to check
#' @param valid_types any of these types are okay.
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd


#' Check inputs for orsf_control_cph()
#'
#' @inheritParams orsf_control_cph
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
#'


#' Check inputs for orsf_control_net()
#'
#' @inheritParams orsf_control_net
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
#'

#' Check inputs for orsf()
#'
#' @inheritParams orsf
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
#'


#' Check inputs for orsf_pd()
#'
#' @inheritParams orsf_pd
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
#'


#' New data have same names as reference data
#'
#' @param new_data data.frame to check
#' @param ref_names character vector of names from reference data
#' @param label_new what to call the new data if an error is printed
#' @param label_ref what to call the reference data if an error is printed.
#' @param check_new_in_ref T/F; make sure all new names are in reference data?
#' @param check_ref_in_new T/F; make sure all reference names are in new data?
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd

check_new_data_names <- function(new_data,
                                 ref_names,
                                 label_new,
                                 label_ref,
                                 check_new_in_ref = FALSE,
                                 check_ref_in_new = TRUE){

 new_names <- names(new_data)

 list_new <- FALSE

 if(check_new_in_ref) list_new <- !(new_names %in% ref_names)

 list_ref <- FALSE

 if(check_ref_in_new) list_ref <- !(ref_names %in% new_names)

 error_new <- any(list_new)
 error_ref <- any(list_ref)

 if(error_new){
  out_msg_new <- paste(
   label_new, " have columns not contained in ", label_ref, ": ",
   paste_collapse(new_names[list_new], last = ' and ')
  )
 }

 if(error_ref){
  out_msg_ref <- paste(
   label_ref, " have columns not contained in ", label_new, ": ",
   paste_collapse(ref_names[list_ref], last = ' and ')
  )
 }

 if(error_new && error_ref){
  out_msg <- c(out_msg_new, '\n Also, ', out_msg_ref)
 }

 if (error_new && !error_ref) {
  out_msg <- c(out_msg_new)
 }

 if (!error_new && error_ref){
  out_msg <- c(out_msg_ref)
 }

 any_error <- error_new | error_ref

 if(any_error){
  stop(out_msg, call. = FALSE)
 }

}

#' New data have same types as reference data
#'
#' If new data have an integer vector where the ref data
#'  had a factor vector, orsf_predict() will yell! Also
#'  it is good practice to make sure users are supplying
#'  consistent data types.
#'
#' @inheritParams check_new_data_names
#' @param ref_types the types of variables in reference data
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
check_new_data_types <- function(new_data,
                                 ref_names,
                                 ref_types,
                                 label_new,
                                 label_ref){

 var_types <- vector(mode = 'character', length = length(ref_names))

 for(i in seq_along(ref_names)){
  var_types[i] <- class(new_data[[ ref_names[i] ]])[1]
 }

 bad_types <- which(var_types != ref_types)

 if(!is_empty(bad_types)){

  vars_to_list <- ref_names[bad_types]
  types_to_list <- var_types[bad_types]

  meat <- paste0('<', vars_to_list, '> has type <',
                 types_to_list, '>', " in ", label_new,
                 "; type <", ref_types[bad_types], "> in ",
                 label_ref, collapse = '\n')

  msg <- paste("some variables in", label_new,
               "have different type in",
               label_ref, ":\n", meat)

  stop(msg, call. = FALSE)

 }

}

#' Check factor variables in new data
#'
#' Factors may have new levels in the testing data, which
#'   would certainly mess up orsf_predict(). So ask user
#'   to fix the factor level before calling predict.
#'
#' @param new_data data that are used for predicting risk or survival
#' @param names_x names of the x variables in training data
#' @param fi_ref factor info from training data
#' @param label_new what to call the new_data if error message is printed.
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd

check_new_data_fctrs <- function(new_data,
                                 names_x,
                                 fi_ref,
                                 label_new){

 fctr_check(new_data, names_x)

 fi_new <- fctr_info(new_data, names_x)

 for(fi_col in fi_ref$cols){
  fctr_check_levels(ref = fi_ref$lvls[[fi_col]],
                    new = fi_new$lvls[[fi_col]],
                    name = fi_col,
                    label_ref = "training data",
                    label_new = label_new)
 }

}





#' check units
#'
#' @param new_data new data to check units in
#'
#' @param ui_train unit information in training data
#'
#' @return nada
#'
#' @noRd
# nocov start
check_units <- function(new_data, ui_train) {

 ui_new <- unit_info(data = new_data, .names = names(ui_train))

 ui_missing <- setdiff(names(ui_train), names(ui_new))

 if(!is_empty(ui_missing)){

  if(length(ui_missing) == 1){
   stop(ui_missing, " had unit attributes in training data but",
        " did not have unit attributes in testing data.",
        " Please ensure that variables in new data have the same",
        " units as their counterparts in the training data.",
        call. = FALSE)
  }

  stop(length(ui_missing),
       " variables (",
       paste_collapse(ui_missing, last = ' and '),
       ") had unit attributes in training",
       " data but did not have unit attributes in new data.",
       " Please ensure that variables in new data have the same",
       " units as their counterparts in the training data.",
       call. = FALSE)

 }

 for(i in names(ui_train)){

  if(ui_train[[i]]$label != ui_new[[i]]$label){

   msg <- paste("variable", i, 'has unit', ui_train[[i]]$label,
                'in the training data but has unit', ui_new[[i]]$label,
                'in new data')

   stop(msg, call. = FALSE)

  }

 }

}
# nocov end

#' Run prediction checks
#'
#' The intent of this function is to protect users from common
#'   inconsistencies that can occur between training data and
#'   testing data. Factor levels in training data need to be
#'   present in the corresponding factors of the testing data,
#'   and all variables used by the model need to be present in
#'   the testing data as well. In addition, the inputs of the
#'   orsf_predict() function are checked.
#'
#' @inheritParams orsf_predict
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd



check_oobag_fun <- function(oobag_fun){

 oobag_fun_args <- names(formals(oobag_fun))

 if(length(oobag_fun_args) != 3) stop(
  "oobag_fun should have 3 input arguments but instead has ",
  length(oobag_fun_args),
  call. = FALSE
 )

 if(oobag_fun_args[1] != 'y_mat') stop(
  "the first input argument of oobag_fun should be named 'y_mat' ",
  "but is instead named '", oobag_fun_args[1], "'",
  call. = FALSE
 )

 if(oobag_fun_args[2] != 'w_vec') stop(
  "the second input argument of oobag_fun should be named 'w_vec' ",
  "but is instead named '", oobag_fun_args[1], "'",
  call. = FALSE
 )

 if(oobag_fun_args[3] != 's_vec') stop(
  "the third input argument of oobag_fun should be named 's_vec' ",
  "but is instead named '", oobag_fun_args[2], "'",
  call. = FALSE
 )

 test_time <- seq(from = 1, to = 5, length.out = 100)
 test_status <- rep(c(0,1), each = 50)

 .y_mat <- cbind(time = test_time, status = test_status)
 .w_vec <- rep(1, times = 100)
 .s_vec <- seq(0.9, 0.1, length.out = 100)

 test_output <- try(oobag_fun(y_mat = .y_mat,
                              w_vec = .w_vec,
                              s_vec = .s_vec),
                    silent = FALSE)

 if(is_error(test_output)){

  stop("oobag_fun encountered an error when it was tested. ",
       "Please make sure your oobag_fun works for this case:\n\n",
       "test_time <- seq(from = 1, to = 5, length.out = 100)\n",
       "test_status <- rep(c(0,1), each = 50)\n\n",
       "y_mat <- cbind(time = test_time, status = test_status)\n",
       "w_vec <- rep(1, times = 100)\n",
       "s_vec <- seq(0.9, 0.1, length.out = 100)\n\n",
       "test_output <- oobag_fun(y_mat = y_mat, w_vec = w_vec, s_vec = s_vec)\n\n",
       "test_output should be a numeric value of length 1",
       call. = FALSE)

 }

 if(!is.numeric(test_output)) stop(
  "oobag_fun should return a numeric output but instead returns ",
  "output of type ", class(test_output)[1],
  call. = FALSE
 )

 if(length(test_output) != 1) stop(
  "oobag_fun should return output of length 1 instead returns ",
  "output of length ", length(test_output),
  call. = FALSE
 )

}

check_beta_fun <- function(beta_fun){

 beta_fun_args <- names(formals(beta_fun))

 if(length(beta_fun_args) != 3) stop(
  "beta_fun should have 3 input arguments but instead has ",
  length(beta_fun_args),
  call. = FALSE
 )

 arg_names_expected <- c("x_node",
                         "y_node",
                         "w_node")

 arg_names_refer <- c('first', 'second', 'third')

 for(i in seq_along(arg_names_expected)){
  if(beta_fun_args[i] != arg_names_expected[i])
   stop(
    "the ", arg_names_refer[i], " input argument of beta_fun ",
    "should be named '", arg_names_expected[i],"' ",
    "but is instead named '", beta_fun_args[i], "'",
    call. = FALSE
   )
 }

 .x_node <- matrix(seq(-1, 1, length.out = 300), ncol = 3)

 test_time <- seq(from = 1, to = 5, length.out = 100)
 test_status <- rep(c(0,1), each = 50)
 .y_node <- cbind(time = test_time, status = test_status)

 .w_node <- matrix(rep(c(1,2,3,4), each = 25), ncol = 1)

 test_output <- try(beta_fun(.x_node, .y_node, .w_node),
                    silent = FALSE)

 if(is_error(test_output)){

  stop("beta_fun encountered an error when it was tested. ",
       "Please make sure your beta_fun works for this case:\n\n",
       ".x_node <- matrix(seq(-1, 1, length.out = 300), ncol = 3)\n\n",
       "test_time <- seq(from = 1, to = 5, length.out = 100)\n",
       "test_status <- rep(c(0,1), each = 50)\n",
       ".y_node <- cbind(time = test_time, status = test_status)\n\n",
       ".w_node <- matrix(rep(c(1,2,3,4), each = 25), ncol = 1)\n\n",
       "test_output <- beta_fun(.x_node, .y_node, .w_node)\n\n",
       "test_output should be a numeric matrix with 1 column and",
       " with nrow(test_output) = ncol(.x_node)",
       call. = FALSE)

 }

 if(!is.matrix(test_output)) stop(
  "beta_fun should return a matrix output but instead returns ",
  "output of type ", class(test_output)[1],
  call. = FALSE
 )

 if(ncol(test_output) != 1) stop(
  "beta_fun should return a matrix with 1 column but instead ",
  " returns a matrix with ", ncol(test_output), " columns.",
  call. = FALSE
 )

 if(nrow(test_output) != ncol(.x_node)) stop(
  "beta_fun should return a matrix with 1 row for each column in x_node ",
  "but instead returns a matrix with ", nrow(test_output), " rows ",
  "in a testing case where x_node has ", ncol(.x_node), " columns",
  call. = FALSE
 )

}

#' check complete cases in new data
#'
#' @param cc (_integer vector_) the indices of complete cases
#' @param na_action the action to be taken for missing values
#'   (see orsf_predict)
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#' @noRd
#'

check_complete_cases <- function(cc, na_action, n_total){

 if(length(cc) != n_total && na_action == 'fail'){
  stop("Please remove missing values from new_data, or impute them.",
       call. = FALSE)
 }

 if(length(cc) == 0){
  stop("There are no observations in new_data with complete data ",
       "for the predictors used by this orsf object.",
       call. = FALSE)
 }

}

