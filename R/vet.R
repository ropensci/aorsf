

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
vet_y <- function(y){

 # check time

 check_arg_type(arg_value = y[, 1],
                arg_name = "time to event",
                expected_type = 'numeric')

 check_arg_gt(arg_value = y[, 1],
              arg_name = "time to event",
              bound = 0)

 # check status

 if(is.logical(y[, 2])) y[, 2] <- as.numeric(y[, 2])

 check_arg_type(arg_value = y[, 2],
                arg_name = "status indicator",
                expected_type = 'numeric')

 status_uni <- unique(y[, 2])

 # a curious case where nothing is censored
 if(all(status_uni == 1)) return(y)

 # status values are modified if they are not all 0 and 1
 if(!is_equivalent(c(0, 1), status_uni)){

  # assume the lowest value of status indicates censoring
  censor_indicator <- min(status_uni)

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

  other_events <- which(y[, 2] > event_indicator)

  # we don't handle other events yet, consider them censored.
  if(!is_empty(other_events)){
   y[other_events, 2] <- censor_indicator
  }

  # The status column needs to contain only 0s and 1s when it
  # gets passed into the C++ routines.
  y[, 2] <- y[, 2] - censor_indicator

  # check again to make sure we don't send something into C++
  # that is going to make the R session crash.
  status_uni <- unique(y[, 2])

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

 y

}
