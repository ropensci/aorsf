

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

 UseMethod('select_cols')

}

# standard data.frame column selection,
# drop = FALSE for type consistency, i.e., always return a data.frame
select_cols.data.frame <- function(data, col_names){

 data[, col_names, drop = FALSE]

}

# data.table doesn't allow column selection in the same manner as data.frames
select_cols.data.table <- function(data, col_names){

 data[, .SD, .SDcols = col_names]

}
