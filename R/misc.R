


is_empty <- function(x) length(x) == 0

is_error <- function(x) inherits(x, 'try-error')

is_cat <- function(x) is.factor(x) & !is.ordered(x)

list_init <- function(.names){
 out <- vector(mode = 'list', length = length(.names))
 names(out) <- .names
 out
}


last_value <- function(x) x[length(x)]

is_aorsf <- function(object){
 inherits(object, 'aorsf')
}
