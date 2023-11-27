
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
check_arg_is_valid <- function(arg_value, arg_name, valid_options,
                               context = NULL) {

 valid_arg <- arg_value %in% valid_options

 if (!valid_arg) {

  expected_values <- paste_collapse(x = valid_options,
                                    sep = ', ',
                                    last = ' or ')

  arg_values <- paste_collapse(x = arg_value,
                               sep = ', ',
                               last = ' or ')

  # context needs an extra bit of text if it isn't null.
  if(!is.null(context)) context <- paste0(" for ", context)

  error_msg <- paste0(
   arg_name, " should be <", expected_values, ">", context,
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


