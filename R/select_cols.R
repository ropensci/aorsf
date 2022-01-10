

select_cols <- function(data, col_names){

 UseMethod('select_cols')

}

select_cols.data.frame <- function(data, col_names){

 data[, col_names, drop = FALSE]

}

# data.table doesn't allow column selection in the same manner as data.frames
select_cols.data.table <- function(data, col_names){

 data[, .SD, .SDcols = col_names]

}
