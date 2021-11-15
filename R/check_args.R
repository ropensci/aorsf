



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

 if(is.null(expected_uni)) return(invisible())

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
   paste0(arg_name, "should contain values of ", expected_values,
          " but has values of ", invalid_values)

  stop(error_msg, call. = FALSE)


 }

}

check_arg_length <- function(arg_value, arg_name, expected_length){

 if(is.null(expected_length)) return(invisible())

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

check_arg_bounds <- function(arg_value, arg_name, bound_lwr, bound_upr){

 arg_value <- arg_value[!is.na(arg_value)]
 if(!is.null(bound_lwr)) check_bound_lwr(arg_value, arg_name, bound_lwr)
 if(!is.null(bound_upr)) check_bound_upr(arg_value, arg_name, bound_upr)

}

check_arg_gt <- function(arg_value, arg_name, bound) {

 if(any(arg_value <= bound)){

  if(length(arg_value) == 1){

   error_msg <-
    paste0(arg_name, " = ", arg_value, "should be >= ", bound)

  } else {

   first_offense <- min(which(arg_value <= bound))

   error_msg <- paste0(arg_name, " should be > ", bound, " but has",
                       " at least one value that is <= ", bound,
                       " (see ", arg_name, "[", first_offense, "])")
  }

  stop(error_msg, call. = FALSE)

 }

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

check_bound_lwr <- function(arg_value, arg_name, bound_lwr) {

 if(any(arg_value < bound_lwr)){

  if(length(arg_value) == 1){

   error_msg <-
    paste0(arg_name, " = ", arg_value, "should be >= ", bound_lwr)

  } else {

   first_offense <- min(which(arg_value < bound_lwr))

   error_msg <- paste0(arg_name, " should be >= ", bound_lwr, " but has",
                       " at least one value that is < ", bound_lwr,
                       " (see ", arg_name, "[", first_offense, "])")
  }

  stop(error_msg, call. = FALSE)

 }

}

check_bound_upr <- function(arg_value, arg_name, bound_upr) {

 if(any(arg_value > bound_upr)){

  if(length(arg_value) == 1){

   error_msg <-
    paste0(arg_name, " = ", arg_value, " should be <= ", bound_upr)

  } else {

   first_offense <- min(which(arg_value > bound_upr))

   error_msg <- paste0(arg_name, " should be <= ", bound_upr, " but has",
                       " at least one value that is > ", bound_upr,
                       " (see ", arg_name, "[", first_offense, "])")
  }

  stop(error_msg, call. = FALSE)
 }

}

check_call <- function(call, expected){

 arg_names <- setdiff( names(call), '' )

 #browser()
 n_frames <- length(sys.frames())

 for (arg_name in arg_names ){

  object_found <- FALSE

  n <- 1

  while(n <= n_frames & !object_found){

   arg_value <- try(
    eval(call[[arg_name]], envir = parent.frame(n = n)),
    silent = TRUE
   )

   if(!inherits(arg_value, 'try-error')) object_found <- TRUE

   n <- n + 1

  }

  if(inherits(arg_value, 'try-error'))
   stop("object '", deparse(call[[arg_name]]),"' not found",
        call. = FALSE)

  if(is.null(arg_value)) return(invisible())

  expected_type <- expected[[arg_name]]$type
  expected_integer <- expected[[arg_name]]$integer
  expected_length <- expected[[arg_name]]$length
  expected_uni <- expected[[arg_name]]$uni
  bound_lwr = expected[[arg_name]]$lwr
  bound_upr = expected[[arg_name]]$upr
  expected_options = expected[[arg_name]]$options
  expected_class = expected[[arg_name]]$class

  if(!is.null(expected_type))
   check_arg_type(arg_name = arg_name,
                  arg_value = arg_value,
                  expected_type = expected_type)

  if(!is.null(expected_integer))
   check_arg_is_integer(arg_name = arg_name,
                        arg_value = arg_value)

  if(!is.null(expected_length))
   check_arg_length(arg_name = arg_name,
                    arg_value = arg_value,
                    expected_length = expected_length)

  if(!is.null(expected_uni))
   check_arg_uni(arg_name = arg_name,
                 arg_value = arg_value,
                 expected_uni = expected_uni)

  if(!is.null(bound_lwr) | !is.null(bound_upr))
   check_arg_bounds(arg_name = arg_name,
                    arg_value = arg_value,
                    bound_lwr = bound_lwr,
                    bound_upr = bound_upr)

  if(!is.null(expected_options))
   check_arg_is_valid(arg_name = arg_name,
                      arg_value = arg_value,
                      valid_options = expected_options)

  if(!is.null(expected_class))
   check_arg_is(arg_name = arg_name,
                arg_value = arg_value,
                expected_class = expected_class)

 }

}


check_var_types <- function(data, valid_types){

 var_types <- purrr::map_chr(data, ~ class(.x)[1])

 good_vars <- var_types %in% valid_types

 if(!all(good_vars)){

  bad_vars <- which(!good_vars)

  vars_to_list <- names(var_types)[bad_vars]
  types_to_list <- var_types[vars_to_list]

  meat <- paste0('<', vars_to_list, '> has type <',
                 types_to_list, '>', collapse = '\n')

  msg <- paste("some variables have unsupported type:\n",
               meat, '\n supported types are', list_things(valid_types)
  )

  stop(msg, call. = FALSE)

 }

}


