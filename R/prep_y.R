

#' Vet the outcome variable
#'
#' Coerce the outcome to be compatible with C++ routines.
#' If there is no feasible way to make it work, throw an error.
#' If there aren't any problems, return the possibly modified outcome.
#'
#' @param y the outcome. For survival, y is a matrix with two columns:
#'   - first column: time values
#'   - second column: status values
#'
#' @return the outcome, possibly modified
#'
#' @noRd
#'
prep_y_surv <- function(data, cols, run_checks = TRUE){

 y <- select_cols(data, cols)

 # if a surv object was included in the formula, it probably
 # doesn't need to be checked here, but it's still checked
 # just to be safe b/c if there is a problem with y it is
 # likely that orsf_fit() will cause R to crash.
 if(length(cols) == 1 && inherits(y[[1]], 'Surv')){
  y <- as.data.frame(as.matrix(y))
  cols <- names(y)
 }

 for(i in cols){
  if(has_units(y[[i]])) y[[i]] <- as.numeric(y[[i]])
 }

 # check time

 if(run_checks){

  check_arg_type(arg_value = y[[1]],
                 arg_name = "time to event",
                 expected_type = 'numeric')

  check_arg_gt(arg_value = y[[1]],
               arg_name = "time to event",
               bound = 0)

  # check status
  if(is.factor(y[[2]]) || is.logical(y[[2]])){
   y[[2]] <- as.numeric(y[[2]])
  }

  check_arg_type(arg_value = y[[2]],
                 arg_name = "status indicator",
                 expected_type = 'numeric')

 }

 # always run so that y mats with 1/2 get passed to CPP as 0/1

 status_uni <- collapse::funique(y[[2]])

 # if nothing is censored
 if(all(status_uni == 1)) return(as_matrix(y))

 # status values are modified if they are not all 0 and 1
 if(!is_equivalent(c(0, 1), status_uni)){

  # assume the lowest value of status indicates censoring
  censor_indicator <- collapse::fmin(status_uni)

  if(censor_indicator != as.integer(censor_indicator)){
   stop("the censoring indicator is not integer valued. ",
        "This can occur when the time column is swapped ",
        "with the status column by mistake. ",
        "Did you enter `status + time` when you meant ",
        "to put `time + status`?", call. = FALSE)
  }

  # assume the integer value 1 above censor is an event
  event_indicator <- censor_indicator + 1

  if( !(event_indicator %in% status_uni) ){
   stop("there does not appear to be an event indicator ",
        "in the status column. The censoring indicator appears to ",
        "be ", censor_indicator, " but there are no values of ",
        1 + censor_indicator, ", which we would assume to be an ",
        " event indicator. This can occur when the time column is ",
        "swapped with the status column by mistake. Did you enter ",
        "`status + time` when you meant  to enter `time + status`? ",
        "If this problem persists, try setting status to 0 for ",
        "censored observations and 1 for events.", call. = FALSE)
  }

  other_events <- which(y[[2]] > event_indicator)

  # we don't handle competing risks yet
  if(!is_empty(other_events)){

   stop("detected >1 event type in status variable. ",
        "Currently only 1 event type is supported. ",
        call. = FALSE)

   # y[other_events, 2] <- censor_indicator

  }

  # The status column needs to contain only 0s and 1s when it
  # gets passed into the C++ routines.
  y[[2]] <- y[[2]] - censor_indicator

  # check to make sure we don't send something into C++
  # that is going to make the R session crash.
  status_uni <- collapse::funique(y[[2]])

  if(!is_equivalent(c(0, 1), status_uni)){
   stop("could not coerce the status column to values of 0 and 1. ",
        "This can occur when the time column is ",
        "swapped with the status column by mistake. Did you enter ",
        "`status + time` when you meant  to enter `time + status`? ",
        "Please modify your data so that the status values are ",
        "0 and 1, with 0 indicating censored observations and 1 ",
        "indicating events.", call. = FALSE)
  }

 }

 as_matrix(y)

}

prep_y_clsf <- function(data, cols, run_checks = TRUE){

 y <- data[[cols]]

 if(is.factor(y)){
  n_class <- length(levels(y))
  y <- as.numeric(y) - 1
 } else {
  n_class <- length(unique(y))
 }

 expand_y_clsf(as_matrix(y), n_class)

}


prep_y_from_orsf <- function(object){

 switch(

  get_tree_type(object),

  'survival' = prep_y_surv(object$data,
                           get_names_y(object),
                           run_checks = FALSE),

  'classification' = prep_y_clsf(object$data,
                                 get_names_y(object),
                                 run_checks = FALSE)

 )


}
