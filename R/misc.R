

#' Are they the same?
#'
#' @param x first set of things
#' @param y second set of things
#'
#' @return are the two sets the same things? `TRUE` if yes
#'
#' @noRd
is_equivalent <- function(x, y) all(x %in% y) && all(y %in% x)

#' Is it empty?
#'
#' @param x object to check
#'
#' @return `TRUE` if length is 0, `FALSE` otherwise
#'
#' @noRd

is_empty <- function(x) length(x) == 0

#' Did it break?
#'
#' @param x object to check
#'
#' @return `TRUE` if the try command failed, `FALSE` otherwise
#'
#' @noRd

is_error <- function(x) inherits(x, 'try-error')


#' initialize a named list
#'
#' @param .names a character vector of names for the new list.
#'
#' @return a list with length `length(.names)` and with names `names`.
#'   Each item in the list is empty, i.e. NULL.
#'
#' @noRd
#'
#' @example
#' list_init(c("A", "B", "C"))

list_init <- function(.names){
 out <- vector(mode = 'list', length = length(.names))
 names(out) <- .names
 out
}


#' Get the last value of a vector
#'
#' @param x a vector
#'
#' @return the last value of `x`
#'
#' @noRd
#'
#' @example
#' last_value(1:10)

last_value <- function(x){
 x[length(x)]
}

#' Helper function to validate aorsf inputs
#'
#' @param object an object to validate
#'
#' @return TRUE if `object` inherits from the 'aorsf' class, FALSE otherwise.
#'
#' @noRd
#'
#' @example
#' is_aorsf("A")

is_aorsf <- function(object){
 inherits(object, 'ObliqueForest')
}


#' paste lists of things
#'
#' @param x a vector of character values
#' @param sep how to space first length(x) - 1 things
#' @param last how to space the last thing
#'
#' @return a string
#'
#' @noRd
#'
#' @example
#' paste_collapse(x = c("first", "second", "third"), last = ' and ')

paste_collapse <- function(x, sep=', ', last = ' or '){

 if(length(x) == 1) return(paste(x))

 if(length(x) == 2) return(paste0(x[1], last, x[2]))

 paste0(paste(x[-length(x)], collapse = sep), trimws(sep), last, x[length(x)])

}

has_units <- function(x){
 inherits(x, 'units')
}

# Clean up after aorsf is unloaded.
.onUnload <- function (libpath) {
 library.dynam.unload("aorsf", libpath)
}



#' Determine whether object has variable importance estimates
#'
#' @param object an object of class 'ObliqueForest'
#'
#' @return `TRUE` if variable importance estimates are present, `FALSE` otherwise
#'
#' @noRd
#'
contains_vi <- function(object) {!is_empty(object$importance)}



#' beautify time
#'
#' @description
#'  Used to make time printouts more readable with verbose progress.
#'  Based on the beautifyTime function in ranger package.
#'
#' @param seconds time in seconds.
#'
#' @noRd
#'
#' @return a string with formatted times

beautifyTime <- function(seconds) {

 result <- ""

 # Add seconds, minutes, hours, days if larger than zero
 out_seconds <- seconds %% 60
 result <- paste(out_seconds, "seconds")

 out_minutes <- (seconds %/% 60) %% 60
 if (seconds %/% 60 == 0) {
  return(result)
 } else if (out_minutes == 1) {
  result <- paste("1 minute,", result)
 } else {
  result <- paste(out_minutes, "minutes,", result)
 }

 out_hours <- (seconds %/% 3600) %% 24
 if (seconds %/% 3600 == 0) {
  return(result)
 } else if (out_hours == 1) {
  result <- paste("1 hour,", result)
 } else {
  result <- paste(out_hours, "hours,", result)
 }

 out_days <- seconds %/% 86400
 if (out_days == 0) {
  return(result)
 } else if (out_days == 1) {
  result <- paste("1 day,", result)
 } else {
  result <- paste(out_days, "days,", result)
 }

 return(result)
}
