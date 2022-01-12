



paste_collapse <- function(x, sep=', ', last = ' or '){

 if(length(x) == 1) return(paste(x))

 if(length(x) == 2) return(paste0(x[1], last, x[2]))

 paste0(paste(x[-length(x)], collapse = sep), trimws(sep), last, x[length(x)])

}


check_arg_type <- function(arg_value, arg_name, expected_type){

 if('numeric' %in% expected_type)
  expected_type <- c(setdiff(expected_type, 'numeric'),
                     'double', 'integer')

 #if(expected_type == 'numeric') expected_type <- c('double', 'integer')

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

check_arg_uni <- function(arg_value, arg_name, expected_uni){

 uni <- unique(arg_value)

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

check_arg_gt <- function(arg_value, arg_name, bound, append_to_msg = NULL){
 check_arg_bound(arg_value,
                 arg_name,
                 bound,
                 relational_operator = 'gt',
                 append_to_msg)
}

check_arg_lt <- function(arg_value, arg_name, bound, append_to_msg = NULL){
 check_arg_bound(arg_value,
                 arg_name,
                 bound,
                 relational_operator = 'lt',
                 append_to_msg)
}

check_arg_gteq <- function(arg_value, arg_name, bound, append_to_msg = NULL){
 check_arg_bound(arg_value,
                 arg_name,
                 bound,
                 relational_operator = 'gteq',
                 append_to_msg)
}

check_arg_lteq <- function(arg_value, arg_name, bound, append_to_msg = NULL){
 check_arg_bound(arg_value,
                 arg_name,
                 bound,
                 relational_operator = 'lteq',
                 append_to_msg)
}

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

check_var_types <- function(data, .names, valid_types){

 var_types <- vector(mode = 'character', length = length(.names))

 for(i in seq_along(.names)){
  var_types[i] <- class(data[[ .names[i] ]])[1]
 }

 good_vars <- var_types %in% valid_types

 if(!all(good_vars)){

  bad_vars <- which(!good_vars)

  vars_to_list <- .names[bad_vars]
  types_to_list <- var_types[bad_vars]

  meat <- paste0('<', vars_to_list, '> has type <',
                 types_to_list, '>', collapse = '\n')

  msg <- paste("some variables have unsupported type:\n",
               meat, '\n supported types are',
               paste_collapse(valid_types, last = ' and '))

  stop(msg, call. = FALSE)

 }

 var_types

}


