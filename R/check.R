



#' @srrstats {G5.2a} *messages produced here (e.g., with `stop()`, `warning()`, `message()`) are unique and make effort to highlight the specific data elements that cause the error*

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


#' strict checks for inputs
#'
#' @srrstats {G2.0} *The function check_arg_length is used to vet inputs.*
#' @srrstats {G2.0a} *Error messages indicate expected lengths*
#'
#' @srrstats {G2.2} *prohibits submission of multivariate input to parameters expected to be univariate.*
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
#' @srrstats {G2.3a} *Using %in% instead of match.arg to allow customized error message when inputs are invalid.*
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
#' @srrstats {G2.4a}  *My design philosophy is that the user should be made aware if their input is the wrong type and the user should also be responsible for correcting the input. Therefore, this function communicates exactly what the user should do to fix the issue.*
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


#' Check variable types
#'
#' @srrstats {G2.1} *types of variables in input data are vetted here.*
#'
#' @srrstats {G2.4} * I think the user should be  aware of type inconsistencies in their data and the user should  also take full responsibility for managing variable types. Therefore, this function communicates problems with variable types to the user but does not fix those problems for the user. It does tell the user exactly how to fix them.*
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
check_control_cph <- function(method, eps, iter_max, pval_max, do_scale){

 check_arg_type(arg_value = method,
                arg_name = 'method',
                expected_type = 'character')

 check_arg_is_valid(arg_value = method,
                    arg_name = 'method',
                    valid_options = c("breslow", "efron"))


 check_arg_type(arg_value = eps,
                arg_name = 'eps',
                expected_type = 'numeric')

 check_arg_gt(arg_value = eps,
              arg_name = 'eps',
              bound = 0)

 check_arg_length(arg_value = eps,
                  arg_name = 'eps',
                  expected_length = 1)


 check_arg_type(arg_value = iter_max,
                arg_name = 'iter_max',
                expected_type = 'numeric')

 check_arg_is_integer(arg_value = iter_max,
                      arg_name = 'iter_max')

 check_arg_gteq(arg_value = iter_max,
                arg_name = 'iter_max',
                bound = 1)

 check_arg_length(arg_value = iter_max,
                  arg_name = 'iter_max',
                  expected_length = 1)


 check_arg_type(arg_value = pval_max,
                arg_name = 'pval_max',
                expected_type = 'numeric')

 check_arg_gt(arg_value = pval_max,
              arg_name = 'pval_max',
              bound = 0)

 check_arg_lteq(arg_value = pval_max,
                arg_name = 'pval_max',
                bound = 1)

 check_arg_length(arg_value = pval_max,
                  arg_name = 'pval_max',
                  expected_length = 1)


 check_arg_type(arg_value = do_scale,
                arg_name = 'do_scale',
                expected_type = 'logical')

 check_arg_length(arg_value = do_scale,
                  arg_name = 'do_scale',
                  expected_length = 1)

 if(!do_scale && iter_max > 1){
  stop("do_scale must be TRUE when iter_max > 1",
       call. = FALSE)
 }

}


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
check_control_net <- function(alpha, df_target){

 check_arg_type(arg_value = alpha,
                arg_name = 'alpha',
                expected_type = 'numeric')

 check_arg_gteq(arg_value = alpha,
                arg_name = 'alpha',
                bound = 0)

 check_arg_lteq(arg_value = alpha,
                arg_name = 'alpha',
                bound = 1)

 check_arg_length(arg_value = alpha,
                  arg_name = 'alpha',
                  expected_length = 1)

 if(!is.null(df_target)){

  check_arg_type(arg_value = df_target,
                 arg_name = 'df_target',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = df_target,
                       arg_name = 'df_target')

 }

}

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
check_orsf_inputs <- function(data_train,
                              formula,
                              control,
                              n_tree,
                              n_split,
                              n_retry,
                              mtry,
                              leaf_min_events,
                              leaf_min_obs,
                              split_min_events,
                              split_min_obs,
                              oobag_pred,
                              oobag_time,
                              oobag_eval_every,
                              importance,
                              tree_seeds,
                              attach_data){

 if(!is.null(data_train)){

  check_arg_is(arg_value = data_train,
               arg_name = 'data_train',
               expected_class = 'data.frame')

  # Minimum event numbers are checked later.
  # Also, later we check to make sure there are at least 2 columns.
  # We specify ncol > 0 here to make the error message that users will
  # receive more specific.
  if(nrow(data_train) == 0 || ncol(data_train) ==  0){
   stop("training data are empty",
        call. = FALSE)
  }

 }

 #' @srrstats {G2.9} issue diagnostic messages for blank column names. In this case, no fixes are applied. Instead, the user is notified by an error message.

 # check for blanks first b/c the check for non-standard symbols
 # will detect blanks with >1 empty characters

 blank_names <- grepl(pattern = '^\\s*$',
                      x = names(data_train))

 if(any(blank_names)){

  s_if_plural_blank_otherwise <- ""

  to_list <- which(blank_names)

  if(length(to_list) > 1) s_if_plural_blank_otherwise <- "s"

  last <- ifelse(length(to_list) == 2, ' and ', ', and ')

  stop("Blank or empty names detected in training data: see column",
       s_if_plural_blank_otherwise, " ",
       paste_collapse(x = to_list, last = last),
       call. = FALSE)

 }

 ns_names <- grepl(pattern = '[^a-zA-Z0-9\\.\\_]+',
                   x = names(data_train))

 if(any(ns_names)){

  last <- ifelse(sum(ns_names) == 2, ' and ', ', and ')

  stop("Non-standard names detected in training data: ",
       paste_collapse(x = names(data_train)[ns_names],
                      last = last),
       call. = FALSE)

 }


 if(!is.null(formula)){

  check_arg_is(arg_value = formula,
               arg_name = 'formula',
               expected_class = 'formula')

  if(length(formula) != 3){
   stop("formula must be two sided, i.e. left side ~ right side",
        call. = FALSE)
  }

  formula_deparsed <- deparse(formula[[3]])

  for( symbol in c("*", "^", ":", "(", ")", "["," ]", "|", "%") ){

   if(grepl(symbol, formula_deparsed, fixed = TRUE)){

    stop("unrecognized symbol in formula: ", symbol,
         "\norsf recognizes '+', '-', and '.' symbols.",
         call. = FALSE)

   }

  }

 }

 check_arg_is(arg_value = control,
              arg_name = 'control',
              expected_class = 'aorsf_control')

 if(!is.null(n_tree)){

  check_arg_type(arg_value = n_tree,
                 arg_name = 'n_tree',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = n_tree,
                       arg_name = 'n_tree')

  check_arg_gteq(arg_value = n_tree,
                 arg_name = 'n_tree',
                 bound = 1)

  check_arg_length(arg_value = n_tree,
                   arg_name = 'n_tree',
                   expected_length = 1)

 }

 if(!is.null(n_split)){

  check_arg_type(arg_value = n_split,
                 arg_name = 'n_split',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = n_split,
                       arg_name = 'n_split')

  check_arg_gteq(arg_value = n_split,
                 arg_name = 'n_split',
                 bound = 1)

  check_arg_length(arg_value = n_split,
                   arg_name = 'n_split',
                   expected_length = 1)

 }

 if(!is.null(n_retry)){

  check_arg_type(arg_value = n_retry,
                 arg_name = 'n_retry',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = n_retry,
                       arg_name = 'n_retry')

  check_arg_gteq(arg_value = n_retry,
                 arg_name = 'n_retry',
                 bound = 0)

  check_arg_length(arg_value = n_retry,
                   arg_name = 'n_retry',
                   expected_length = 1)

 }

 if(!is.null(mtry)){

  check_arg_type(arg_value = mtry,
                 arg_name = 'mtry',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_name = 'mtry',
                       arg_value = mtry)

  check_arg_gteq(arg_name = 'mtry',
                 arg_value = mtry,
                 bound = 2)

  check_arg_length(arg_name = 'mtry',
                   arg_value = mtry,
                   expected_length = 1)

 }

 if(!is.null(leaf_min_events)){

  check_arg_type(arg_value = leaf_min_events,
                 arg_name = 'leaf_min_events',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = leaf_min_events,
                       arg_name = 'leaf_min_events')

  check_arg_gteq(arg_value = leaf_min_events,
                 arg_name = 'leaf_min_events',
                 bound = 1)

  check_arg_length(arg_value = leaf_min_events,
                   arg_name = 'leaf_min_events',
                   expected_length = 1)
 }

 if(!is.null(leaf_min_obs)){

  check_arg_type(arg_value = leaf_min_obs,
                 arg_name = 'leaf_min_obs',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = leaf_min_obs,
                       arg_name = 'leaf_min_obs')

  check_arg_gteq(arg_value = leaf_min_obs,
                 arg_name = 'leaf_min_obs',
                 bound = 1)

  check_arg_length(arg_value = leaf_min_obs,
                   arg_name = 'leaf_min_obs',
                   expected_length = 1)

 }

 if(!is.null(split_min_events)){

  check_arg_type(arg_value = split_min_events,
                 arg_name = 'split_min_events',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = split_min_events,
                       arg_name = 'split_min_events')

  check_arg_gteq(arg_value = split_min_events,
                 arg_name = 'split_min_events',
                 bound = 1)

  check_arg_length(arg_value = split_min_events,
                   arg_name = 'split_min_events',
                   expected_length = 1)
 }

 if(!is.null(split_min_obs)){

  check_arg_type(arg_value = split_min_obs,
                 arg_name = 'split_min_obs',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = split_min_obs,
                       arg_name = 'split_min_obs')

  check_arg_gteq(arg_value = split_min_obs,
                 arg_name = 'split_min_obs',
                 bound = 1)

  check_arg_length(arg_value = split_min_obs,
                   arg_name = 'split_min_obs',
                   expected_length = 1)

 }

 if(!is.null(oobag_pred)){

  check_arg_type(arg_value = oobag_pred,
                 arg_name = 'oobag_pred',
                 expected_type = 'logical')

  check_arg_length(arg_value = oobag_pred,
                   arg_name = 'oobag_pred',
                   expected_length = 1)

 }

 if(!is.null(oobag_time)){

  check_arg_type(arg_value = oobag_time,
                 arg_name = 'oobag_time',
                 expected_type = 'numeric')

  check_arg_length(arg_value = oobag_time,
                   arg_name = 'oobag_time',
                   expected_length = 1)

  check_arg_gt(arg_value = oobag_time,
               arg_name = 'oobag_time',
               bound = 0)

 }


 if(!is.null(oobag_eval_every)){

  check_arg_type(arg_value = oobag_eval_every,
                 arg_name = 'oobag_eval_every',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = oobag_eval_every,
                       arg_name = 'oobag_eval_every')

  check_arg_gteq(arg_value = oobag_eval_every,
                 arg_name = 'oobag_eval_every',
                 bound = 1)

  check_arg_lteq(arg_value = oobag_eval_every,
                 arg_name = 'oobag_eval_every',
                 bound = n_tree)

  check_arg_length(arg_value = oobag_eval_every,
                   arg_name = 'oobag_eval_every',
                   expected_length = 1)

 }

 if(!is.null(tree_seeds)){

  check_arg_type(arg_value = tree_seeds,
                 arg_name = 'tree_seed',
                 expected_type = 'numeric')

  check_arg_is_integer(tree_seeds, arg_name = 'tree_seeds')

  if(length(tree_seeds) != n_tree){

   stop('tree_seeds should have length <', n_tree,
        "> (the number of trees) but instead has length <",
        length(tree_seeds), ">", call. = FALSE)

  }

 }

 if(!is.null(attach_data)){

  check_arg_type(arg_value = attach_data,
                 arg_name = 'attach_data',
                 expected_type = 'logical')

  check_arg_length(arg_value = attach_data,
                   arg_name = 'attach_data',
                   expected_length = 1)

 }

}

#' Check inputs for orsf_pd()
#'
#' @inheritParams orsf_pd_summary
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
#'
check_pd_inputs <- function(object,
                            expand_grid = NULL,
                            prob_values = NULL,
                            prob_labels = NULL,
                            oobag = NULL,
                            risk = NULL){

 check_arg_is(arg_value = object,
              arg_name = 'object',
              expected_class = 'aorsf')

 if(!is.null(expand_grid)){

  check_arg_type(arg_value = expand_grid,
                 arg_name = 'expand_grid',
                 expected_type = 'logical')

  check_arg_length(arg_value = expand_grid,
                   arg_name = 'expand_grid',
                   expected_length = 1)

 }

 if(!is.null(prob_values)){

  check_arg_type(arg_value = prob_values,
                 arg_name = 'prob_values',
                 expected_type = 'numeric')

  check_arg_gteq(arg_value = prob_values,
                 arg_name = 'prob_values',
                 bound = 0)

  check_arg_lteq(arg_value = prob_values,
                 arg_name = 'prob_values',
                 bound = 1)

 }

 if(!is.null(prob_labels)){

  check_arg_type(arg_value = prob_labels,
                 arg_name = 'prob_labels',
                 expected_type = 'character')

 }

 if(!is.null(oobag)){

  check_arg_type(arg_value = oobag,
                 arg_name = 'oobag',
                 expected_type = 'logical')

  check_arg_length(arg_value = oobag,
                   arg_name = 'oobag',
                   expected_length = 1)

 }

 if(!is.null(risk)){

  check_arg_type(arg_value = risk,
                 arg_name = 'risk',
                 expected_type = 'logical')

  check_arg_length(arg_value = risk,
                   arg_name = 'risk',
                   expected_length = 1)

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
#' @srrstats {G2.4} *I think the user should be aware of type inconsistencies in their data and the user should  also take full responsibility for managing variable types. Therefore, this function communicates problems with variable types to the user but does not fix those problems for the user. It does tell the user exactly how to fix them though.*
#'
#' @srrstats {G2.6} *Ensure all one-dimensional inputs in new_data have been processed to match their corresponding inputs in the training data for the aorsf model. I.e., if a variable was a factor in the training data, it needs to be a factor with the exact same levels in the testing data. If a variable was an integer in the training data, it must also be an integer in the testing data.*
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

  msg <- paste("some variables in ", label_new,
               " have different type in ",
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


#' check levels of individual factor
#'
#' @param ref levels of factor in reference data
#' @param new levels of factor in new data
#' @param name name of the factor variable
#' @param label_ref what to call reference data if error message is printed.
#' @param label_new what to call new data if error message is printed.
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd

fctr_check_levels <- function(ref,
                              new,
                              name,
                              label_ref,
                              label_new){

 list_new  <- !(new %in% ref)

 if(any(list_new)){

  out_msg <- paste0(
   "variable ", name, " in ", label_new,
   " has levels not contained in ", label_ref, ": ",
   paste_collapse(new[list_new], last = ' and ')
  )

  stop(out_msg, call. = FALSE)

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

check_predict <- function(object, new_data, pred_horizon, risk){

 if(!is.null(new_data)){

  check_arg_is(arg_value = new_data,
               arg_name = 'new_data',
               expected_class = 'data.frame')

  if(nrow(new_data) == 0 || ncol(new_data) ==  0){
   stop("new data are empty",
        call. = FALSE)
  }

 }

 if(!is.null(pred_horizon)){

  check_arg_type(arg_value = pred_horizon,
                 arg_name = 'pred_horizon',
                 expected_type = 'numeric')

  check_arg_gt(arg_value = pred_horizon,
               arg_name = 'pred_horizon',
               bound = 0)

 }

 if(!is.null(risk)){

  check_arg_type(arg_value = risk,
                 arg_name = 'risk',
                 expected_type = 'logical')

  check_arg_length(arg_value = risk,
                   arg_name = 'risk',
                   expected_length = 1)

 }

 ui_train <- get_unit_info(object)

 # check unit info for new data if training data had unit variables
 if(!is_empty(ui_train)) check_units(new_data, ui_train)

 if(any(pred_horizon > get_max_time(object))){

  stop("prediction horizon should ",
       "be <= max follow-up time ",
       "observed in training data: ",
       get_max_time(object),
       call. = FALSE)

 }

 if(!all(order(pred_horizon) == seq(length(pred_horizon)))){
  stop("prediction horizons must be entered in ascending order, e.g.,",
       "pred_horizon = c(5, 10) instead of pred_horizon = c(10, 5)",
       call. = FALSE)
 }

 check_new_data_names(new_data  = new_data,
                      ref_names = get_names_x(object),
                      label_new = "new_data",
                      label_ref = 'training data')

 check_new_data_types(new_data  = new_data,
                      ref_names = get_names_x(object),
                      ref_types = get_types_x(object),
                      label_new = "new_data",
                      label_ref = 'training data')

 check_new_data_fctrs(new_data  = new_data,
                      names_x   = get_names_x(object),
                      fi_ref    = get_fctr_info(object),
                      label_new = "new_data")

 #' @srrstats {G2.6} *ensure that one-dimensional inputs are appropriately pre-processed. aorsf does not deal with missing data as many other R packages are very good at dealing with it.*

 #' @srrstats {G2.13} *check for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*

 #' @srrstats {G2.15} *Never pass data with potential missing values to any base routines.*

 if(any(is.na(new_data[, get_names_x(object)]))){
  stop("Please remove missing values from new_data, or impute them.",
       call. = FALSE)
 }

 #' @srrstats {G2.16} *Throw hard errors if undefined values are detected.*

 for(i in c(get_names_x(object))){

  if(any(is.infinite(new_data[[i]]))){
   stop("Please remove infinite values from ", i, ".",
        call. = FALSE)
  }

  # NaN values trigger is.na(), so this probaly isn't needed.
  # if(any(is.nan(new_data[[i]]))){
  #  stop("Please remove NaN values from ", i, ".",
  #       call. = FALSE)
  # }

 }


}

check_oobag_fun <- function(oobag_fun){

 oobag_fun_args <- names(formals(oobag_fun))

 if(length(oobag_fun_args) != 2) stop(
  "oobag_fun should have 2 input arguments but instead has ",
  length(oobag_fun_args),
  call. = FALSE
 )

 if(oobag_fun_args[1] != 'y_mat') stop(
  "the first input argument of oobag_fun should be named 'y_mat' ",
  "but is instead named '", oobag_fun_args[1], "'",
  call. = FALSE
 )

 if(oobag_fun_args[2] != 's_vec') stop(
  "the second input argument of oobag_fun should be named 's_vec' ",
  "but is instead named '", oobag_fun_args[2], "'",
  call. = FALSE
 )

 test_time <- seq(from = 1, to = 5, length.out = 100)
 test_status <- rep(c(0,1), each = 50)

 .y_mat <- cbind(time = test_time, status = test_status)
 .s_vec <- seq(0.9, 0.1, length.out = 100)

 test_output <- try(oobag_fun(y_mat = .y_mat, s_vec = .s_vec),
                    silent = FALSE)

 if(is_error(test_output)){

  stop("oobag_fun encountered an error when it was tested. ",
       "Please make sure your oobag_fun works for this case:\n\n",
       "test_time <- seq(from = 1, to = 5, length.out = 100)\n",
       "test_status <- rep(c(0,1), each = 50)\n\n",
       "y_mat <- cbind(time = test_time, status = test_status)\n",
       "s_vec <- seq(0.9, 0.1, length.out = 100)\n\n",
       "test_output <- oobag_fun(y_mat = y_mat, s_vec = s_vec)\n\n",
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


