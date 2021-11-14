

fctr_check <- function(data, .names){

 chrs <- c()
 for( .name in .names ){
  if(is.character(data[[.name]])){
   chrs <- c(chrs, .name)
  }
 }

 if(is_empty(chrs)) return(NULL)

 chrs <- paste_collapse(chrs, last = ' and ')

 stop(
  "character variables in data should be converted to factors.\n",
  "Here are the character variables I detected: ", chrs, call.= FALSE
 )

}


fctr_names <- function(data, .names){

 fctrs <- vector(mode = 'character')

 for( .name in .names ){
  if(is.factor(data[[.name]]) & !is.ordered(data[[.name]])){
   fctrs <- c(fctrs, .name)
  }
 }

 fctrs

}


fctr_info <- function(data, .names, fctr_sep = '_'){

 fctr_check(data, .names)

 fctr_info <- vector(mode = 'list', length = 3L)
 names(fctr_info) <- c('cols', 'lvls', 'keys')

 fctr_names <- fctr_names(data, .names)

 # dont waste time if there aren't any factors
 if(is_empty(fctr_names)) return(fctr_info)

 fctr_count <- length(fctr_names)

 fctr_info$cols <- fctr_names
 fctr_info$lvls <- vector(mode = 'list', length = fctr_count)
 fctr_info$keys <- vector(mode = 'list', length = fctr_count)

 names(fctr_info$lvls) <- fctr_names
 names(fctr_info$keys) <- fctr_names

 for( fctr in fctr_names ){

   lvls <- levels(data[[fctr]])
   fctr_info$lvls[[fctr]] <- lvls
   fctr_info$keys[[fctr]] <- paste(fctr, lvls, sep = fctr_sep)


 }

 fctr_info

}

