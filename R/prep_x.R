
prep_x <- function(data,
                   fi,
                   cols,
                   means,
                   standard_deviations){

 cols_numeric <- setdiff(cols, fi$cols)

 x <- ref_code(data, fi, cols)

 for(i in cols_numeric){

  if(has_units(x[[i]])) x[[i]] <- as.numeric(x[[i]])

  # can't modify by reference here, it would modify the user's data
  x[[i]] <- (x[[i]] - means[i]) / standard_deviations[i]

 }

 as_matrix(x)

}


prep_x_from_orsf <- function(object,
                             data = object$data){

 prep_x(data,
        fi = get_fctr_info(object),
        cols = get_names_x(object),
        means = get_means(object),
        standard_deviations = get_standard_deviations(object))



}
