


#' Factor tests
#'
#' @param data data frame to check factors in
#' @param .names names of variables in `data` to consider.
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
#'
fctr_check <- function(data, .names){

 # using select_cols in case data is a data.table object.
 # I take the intersection of .names with names of data b/c
 # I don't want this function to throw an error if the user
 # specified a variable that was not in the dataframe.
 # (there is another check function dedicated to that kind of error)
 chr_index <- which(
  vapply(X = select_cols(data, intersect(.names, names(data))),
         FUN = is.character,
         FUN.VALUE = logical(1),
         USE.NAMES = FALSE)
 )

 if(is_empty(chr_index)) return(NULL)

 chrs <- .names[chr_index]

 chrs <- paste_collapse(chrs, last = ' and ')

 stop(
  "character variables in data should be converted to factors.\n",
  "Here are the character variables I detected: ", chrs, call.= FALSE
 )


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

#' Factor tests
#'
#' There is a good chance that someone using Surv(time, status) ~ .
#'   will forget that inside of the '.' sits an ID variable.
#'   Catching that and sending an informative error will likely
#'   be appreciated.
#'
#' @param data data frame to check factors in
#' @param .names names of variables in `data` to consider.
#'
#' @return check functions 'return' errors and the intent is
#'   to return nothing if nothing is wrong,
#'   so hopefully nothing is returned.
#'
#' @noRd
#'
fctr_id_check <- function(data, .names){

 for(.name in .names) {

  if(is.factor(data[[.name]])){

   if(length(levels(data[[.name]])) == nrow(data)){
    stop("factor variable ", .name, " has as many levels as there ",
         "are rows in the training data. Is ", .name, " an id variable?",
         call. = FALSE)
   }

  }

 }



}


#' Factor information
#'
#' @param data data frame to check factors in
#' @param .names names of variables in `data` to consider.
#' @param fctr_sep how to separate factor variable names from levels.
#'
#' @return a list describing factor variables in `data`.
#'
#' @noRd
#'
#' @example
#' fctr_info(pbc_orsf, .names = c('sex','stage'))

fctr_info <- function(data, .names, fctr_sep = '_'){

 fctr_check(data, .names)

 fctrs <- vector(mode = 'character')
 ordrd <- vector(mode = 'logical')

 for( .name in .names ){
  if(is.factor(data[[.name]])){
   fctrs <- c(fctrs, .name)
   ordrd <- c(ordrd, is.ordered(data[[.name]]))
  }
 }

 fctr_info <- vector(mode = 'list', length = 4L)
 names(fctr_info) <- c('cols', 'lvls', 'keys', 'ordr')

 # dont waste time if there aren't any factors
 if(is_empty(fctrs)) return(fctr_info)

 fctr_info$cols <- fctrs
 fctr_info$ordr <- ordrd
 fctr_info$lvls <- vector(mode = 'list', length = length(fctrs))
 fctr_info$keys <- vector(mode = 'list', length = length(fctrs))

 names(fctr_info$lvls) <- fctrs
 names(fctr_info$keys) <- fctrs

 for( i in seq_along(fctrs) ){

  fctr <- fctrs[i]

  lvls <- levels(data[[fctr]])
  fctr_info$lvls[[fctr]] <- lvls

  if(!fctr_info$ordr[i])
   fctr_info$keys[[fctr]] <- paste(fctr, lvls, sep = fctr_sep)

 }

 fctr_info

}

