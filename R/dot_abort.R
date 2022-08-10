

dot_abort <- function(.dots, .args){

 for(i in seq_along(.dots)){

  .match_indices <- utils::adist(x = .dots[i],
                                 y = .args,
                                 fixed = TRUE,
                                 costs = c(ins = 1,
                                           del = 1,
                                           sub = 2))

  .match_index <- which.min(.match_indices)

  .dots[i] <- paste('  ', .dots[i], ' - did you mean ',
                    .args[.match_index], '?', sep = '')
 }


 stop("the following arguments were not recognized:\n",
      paste(.dots, collapse = '\n'),
      call. = FALSE)

}
