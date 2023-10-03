
# misc functions used for tests

no_miss_list <- function(l){

 sapply(l, function(x){

  if(is.list(x)) {return(no_miss_list(x))}

  any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x))

 })

}

add_noise <- function(x, eps = .Machine$double.eps){

 noise <- rnorm(length(x), mean = 0, sd = eps/2)
 noise <- pmin(noise, eps)
 noise <- pmax(noise, -eps)

 x + noise

}

change_scale <- function(x, mult_by = 1/2){
 x * mult_by
}
