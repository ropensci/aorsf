

fctr_check <- function(data, .names){

 chrs <- c()

 for( .name in .names ){

  if(is.character(data[[.name]])){
   chrs <- c(chrs, .name)
  }

  if(is.factor(data[[.name]])){
   if(length(levels(data[[.name]])) == nrow(data)){
    stop("factor variable ", .name, " has as many levels as there ",
         "are rows in the training data. Is ", .name, "an id variable?")
   }
  }

 }

 if(is_empty(chrs)) return(NULL)

 chrs <- paste_collapse(chrs, last = ' and ')

 stop(
  "character variables in data should be converted to factors.\n",
  "Here are the character variables I detected: ", chrs, call.= FALSE
 )


}


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

