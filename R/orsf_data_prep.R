
orsf_data_prep <- function(data, ...){
 UseMethod('orsf_data_prep')
}

orsf_data_prep.list <- function(data, ...){

 lengths <- vapply(data, length, integer(1))

 if(! all(lengths == lengths[1])){

  length_tbl <- table(lengths)
  length_mode <- as.numeric(names(length_tbl)[which.max(length_tbl)])

  mismatch <- lengths[names(which(lengths != length_mode))]

  mismatch <-
   paste(" -", names(mismatch),
         'has length', mismatch,
         collapse = '\n')

  mismatch <-
   paste(mismatch,
         '\n - all other variables have length ', length_mode,
         sep = '')

  stop("unable to cast data (a list) into a data.frame.\n",
       mismatch, call. = FALSE)

 }

 data <-
  tryCatch(as.data.frame(data), error = function(e) e$message)

 if(!is.data.frame(data)){
  stop("Could not coerce data (a list) into a data.frame object.\n",
       "Running as.data.frame(data) ",
       "produced this error message:\n\"", data, "\"",
       call. = FALSE)
 }

 data

}

orsf_data_prep.recipe <- function(data, ...){

 getElement(data, 'template')

}

orsf_data_prep.data.frame <- function(data, ...){
 data
}
