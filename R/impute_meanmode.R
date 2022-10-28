

#' simple mode function for mean/mode imputation

mode <- function(x){

 # don't use randomness here - it will lead to unexpected differences
 # between orsf and orsf_train when set.seed is used to sync them up
 tab <- table(x)
 modes <- names(tab)[tab == max(tab)][1]
 modes[1]

}

get_means <- function(data, cols){

 # do nothing if there aren't any values in cols
 if(is_empty(cols)) return(NULL)

 lapply(select_cols(data, cols), mean, na.rm = TRUE)

}

get_modes <- function(data, cols){

 # do nothing if there aren't any values in cols
 if(is_empty(cols)) return(NULL)

 lapply(select_cols(data, cols), mode)

}

#' dealing with missing data
impute_meanmode <- function(data, cols, values){
 UseMethod('impute_meanmode')
}

# warning: modifies data in place
impute_meanmode.data.table <- function(data, cols, values){

 for(j in cols){

  if(any(is.na(data[[j]]))){

   set(data,
       i = which(is.na(data[[j]])),
       j = values[[j]])

  }

 }

 data

}

# makes a copy
impute_meanmode.data.frame <- function(data, cols, values){

 for(j in cols){

  if(any(is.na(data[[j]]))){

   data[which(is.na(data[[j]])), j] <- values[[j]]

  }

 }

 data

}
