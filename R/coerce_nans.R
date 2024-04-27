coerce_nans <- function(x, to){
 UseMethod('coerce_nans')
}

coerce_nans.list <- function(x, to){

 lapply(x, coerce_nans, to = to)

}

coerce_nans.array <- coerce_nans.matrix <- function(x, to){

 if(any(is.nan(x))){
  x[is.nan(x)] <- to
 }

 x

}

