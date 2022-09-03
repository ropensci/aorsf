

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
 inherits(object, 'orsf_fit')
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

# Clean up after aorsf is unloaded.
.onUnload <- function (libpath) {
 library.dynam.unload("aorsf", libpath)
}

