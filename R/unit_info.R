


#' get unit info
#'
#' @param data dataset to extract unit info from
#' @param .names a character vector containing names that have unit class.
#'
#' @return a list with unit info
#'
#' @noRd

unit_info <- function(data, .names){

 out <- list()

 if(is_empty(.names)) return(out)

 for(i in .names){

  u <- attr(data[[i]], 'units')

  if(inherits(u, 'symbolic_units')){

   n <- u$numerator
   d <- u$denominator
   l <- paste(n)

   if(!is_empty(d)) l <- paste0(l ,'/', d)

   out[[i]] <- list(numerator = n,
                    denominator = d,
                    label = l)

  }



 }

 out


}
