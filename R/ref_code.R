

#' Reference coding for factors
#'
#' Use reference coding of factors to make sure there is no possibility of
#'   collinearity while I am running coxph routines inside of orsf().
#'
#' @param x_data data frame to convert to one-hot encoding
#' @param fi factor information
#' @param names_x_data names of the variables to consider
#'
#' @return the same data with factors formatted in one-hot encoding.
#'
#' @noRd
#'
#' @examples
#' dat_refcoded <- ref_code(x_data = pbc_orsf,
#'                   fi = fctr_info(pbc_orsf, .names = 'sex'),
#'                   names_x_data = c('age', 'sex'))
#'
#' head(dat_refcoded)
#'
ref_code <- function (x_data, fi, names_x_data){

 # Will use these original names to help re-order the output

 for(i in seq_along(fi$cols)){

  if(fi$cols[i] %in% names(x_data)){

   if(fi$ordr[i]){

    x_data[[ fi$cols[i] ]] <- as.integer( x_data[[ fi$cols[i] ]] )

   } else {

    # make a matrix for each factor
    mat <- matrix(0,
                  nrow = nrow(x_data),
                  ncol = length(fi$lvls[[i]])
    )

    colnames(mat) <- fi$keys[[i]]

    # missing values of the factor become missing rows
    mat[is.na(x_data[[fi$cols[i]]]), ] <- NA_integer_

    # we will one-hot encode the matrix and then bind it to data,
    # replacing the original factor column. Go through the matrix
    # column by column, where each column corresponds to a level
    # of the current factor (indexed by i). Flip the values
    # of the j'th column to 1 whenever the current factor's value
    # is the j'th level.

    for (j in seq(ncol(mat))) {

     # find which rows to turn into 1's. These should be the
     # indices in the currect factor where it's value is equal
     # to the j'th level.
     hot_rows <- which( x_data[[fi$cols[i]]] == fi$lvls[[i]][j] )

     # after finding the rows, flip the values from 0 to 1
     if(!is_empty(hot_rows)){
      mat[hot_rows , j] <- 1
     }

    }

    # data[[fi$cols[i]]] <- NULL

    x_data <- cbind(x_data, mat)

   }

  }



 }

 OH_names <- names_x_data

 for (i in seq_along(fi$cols)){

  if(fi$cols[i] %in% names_x_data){
   if(!fi$ordr[i]){
    OH_names <- insert_vals(
     vec = OH_names,
     where = which(fi$cols[i] == OH_names),
     what = fi$keys[[i]][-1]
    )
   }
  }

 }

 select_cols(x_data, OH_names)

}




#' insert some value(s) into a vector
#'
#'
#' @param vec the vector to be edited
#' @param where where to insert values
#' @param what what is to be inserted
#'
#' @return a vector with new values inserted
#'
#' @details The value that was originally in the vector at
#'   `where` will be replaced (see example).
#'
#' @noRd
#'
#' @examples
#'
#' insert_vals(vec = 1:4, where = 4, what = c(5L))
#'
insert_vals <- function(vec, where, what){

 stopifnot(
  typeof(what) == typeof(vec),
  where >= 1 & where <= length(vec)
 )

 if(where == 1){

  if(length(vec) == 1) return(c(what)) else return(c(what, vec[-1]))
 }

 if(where == length(vec)) return(c(vec[1:(length(vec)-1)], what))

 vec_left <- vec[1:(where-1)]
 vec_right <- vec[(where+1):length(vec)]

 c(vec_left, what, vec_right)

}

