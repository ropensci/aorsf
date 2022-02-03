

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
#' @srrstats {G2.10} *set drop = FALSE to ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behavior, and all column-extraction operations behave consistently regardless of the class of tabular data used as input.*
select_cols.data.frame <- function(data, col_names){

 data[, col_names, drop = FALSE]

}

# data.table doesn't allow column selection in the same manner as data.frames
select_cols.data.table <- function(data, col_names){

 data[, .SD, .SDcols = col_names]

}
