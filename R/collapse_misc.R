
as_matrix <- function(x){
 collapse::qM(x)
}

#' subset data by columns
#'
#' Made this function so that data.table inputs wouldn't require
#'   special treatment
#'
#' @param data data to select columns from
#' @param col_names names of the columns to select
#'
#' @return data[, col_names]
#' @noRd
#'
select_cols <- function(data, col_names){

 collapse::fselect(data, col_names, return = 'data')

}


data_impute <- function(data, cols, values){

 for(col in cols){

  if(anyNA(data[[col]]))
   data <- collapse::replace_NA(data,
                                cols = col,
                                value = values[[col]])

 }

 data

}


data_impute_nocheck <- function(data, cols, values){

 for(col in cols)
  data <- collapse::replace_NA(data,
                               cols = col,
                               value = values[[col]])

 data

}


