

#' Title
#'
#' @param object
#' @param new_data
#' @param times
#'
#' @return
#' @export
#'
#' @examples
predict.aorsf <- function(object, new_data, times, risk = TRUE){

 Call <- match.call()

 check_call(
  Call,
  expected = list(
   'new_data' = list(
    class = 'data.frame'
   ),
   'times' = list(
    type = 'numeric',
    lwr = 0
   ),
   'risk' = list(
    type = 'logical',
    length = 1
   )
  )
 )

 check_new_data_names(new_data  = new_data,
                      ref_names = object$names_x,
                      label_new = deparse(Call$new_data),
                      label_ref = 'training data')

 check_new_data_fctrs(new_data  = new_data,
                      names_x   = object$names_x,
                      fi_ref    = object$fctr_info,
                      label_new = deparse(Call$new_data))

 if(!all(order(times) == seq(length(times)))){
  stop("times must be entered in ascending order, e.g.,",
       "times = c(5, 10) instead of times = c(10, 5)",
       call. = FALSE)
 }

 x_new <- as.matrix(
  one_hot(data = new_data,
          fi = object$fctr_info,
          names_x_data = object$names_x)
 )

 if(length(times) == 1){
  return(orsf_pred_uni(object$forest, x_new, times, risk))
 }

 orsf_pred_multi(object$forest, x_new, times, risk)

}


check_new_data_names <- function(new_data,
                                 ref_names,
                                 label_new,
                                 label_ref,
                                 check_new_in_ref = FALSE,
                                 check_ref_in_new = TRUE){

 new_names <- names(new_data)

 list_new <- FALSE

 if(check_new_in_ref) list_new <- !(new_names %in% ref_names)

 list_ref <- FALSE

 if(check_ref_in_new) list_ref <- !(ref_names %in% new_names)

 error_new <- any(list_new)
 error_ref <- any(list_ref)

 if(error_new){
  out_msg_new <- paste(
   label_new, " have columns not contained in ", label_ref, ": ",
   paste_collapse(new_names[list_new], last = ' and ')
  )
 }

 if(error_ref){
  out_msg_ref <- paste(
   label_ref, " have columns not contained in ", label_new, ": ",
   paste_collapse(ref_names[list_ref], last = ' and ')
  )
 }

 if(error_new && error_ref){
  out_msg <- c(out_msg_new, '\n Also, ', out_msg_ref)
 }

 if (error_new && !error_ref) {
  out_msg <- c(out_msg_new)
 }

 if (!error_new && error_ref){
  out_msg <- c(out_msg_ref)
 }

 any_error <- error_new | error_ref

 if(any_error){
  stop(out_msg, call. = FALSE)
 }

}

check_new_data_fctrs <- function(new_data,
                                 names_x,
                                 fi_ref,
                                 label_new){

 fctr_check(new_data, names_x)

 fi_new <- fctr_info(new_data, names_x)

 for(fi_col in fi_ref$cols){
  fctr_check_levels(ref = fi_ref$lvls[[fi_col]],
                    new = fi_new$lvls[[fi_col]],
                    name = fi_col,
                    label_ref = "training data",
                    label_new = label_new)
 }

}

fctr_check_levels <- function(ref,
                              new,
                              name,
                              label_ref,
                              label_new){

 list_new  <- !(new %in% ref)

 if(any(list_new)){

  out_msg <- paste0(
   "variable ", name, " in ", label_new,
   " has levels not contained in ", label_ref, ": ",
   paste_collapse(new[list_new], last = ' and ')
  )

  stop(out_msg, call. = FALSE)

 }


}

