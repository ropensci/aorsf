



#' print aosrf object
#'
#' @param x an object of class 'aorsf'
#' @param ... not used
#'
#' @return nothing - just print output to console
#' @export
#'
print.aorsf <- function(x, ...){

 info_n_obs           <- get_n_obs(x)
 info_n_events        <- get_n_events(x)
 info_n_tree          <- get_n_tree(x)
 info_max_time        <- get_max_time(x)
 info_n_leaves_mean   <- table.glue::table_value(get_n_leaves_mean(x))
 info_mtry            <- get_mtry(x)
 info_names_x         <- get_names_x(x)
 info_leaf_min_obs    <- get_leaf_min_obs(x)
 info_leaf_min_events <- get_leaf_min_events(x)

 info_oobag_c <- ifelse(test = has_oobag(x),
                        yes = table.glue::table_value(x$eval_oobag$c_harrell),
                        no = "none")


 cat('---------- Oblique random survival forest\n',
     paste0('          N observations: ', info_n_obs          ),
     paste0('                N events: ', info_n_events       ),
     paste0('                 N trees: ', info_n_tree         ),
     paste0('      N predictors total: ', length(info_names_x)),
     paste0('   N predictors per node: ', info_mtry           ),
     paste0(' Average leaves per tree: ', info_n_leaves_mean  ),
     paste0('Min observations in leaf: ', info_leaf_min_obs   ),
     paste0('      Min events in leaf: ', info_leaf_min_events),
     paste0('         OOB C-statistic: ', info_oobag_c        ),
     '\n-----------------------------------------',
     sep = '\n')




}

print.aorsf_summary <- function(x, ...){

 .dots <- list(...)

 if(is.null(.dots$topn))       .dots$topn <- 5
 if(is.null(.dots$nrows))      .dots$nrows <- 100
 if(is.null(.dots$class))      .dots$class <- TRUE
 if(is.null(.dots$row.names))  .dots$row.names <- TRUE
 if(is.null(.dots$col.names))  .dots$col.names <- "auto"
 if(is.null(.dots$print.keys)) .dots$print.keys <- TRUE
 if(is.null(.dots$digits))     .dots$digits <- 3

 .dots$x <- x$data_inputs

 banner_input_length <-
  max(vapply(
   utils::capture.output(do.call(print, .dots)),
   nchar,
   integer(1)
  ))

 banner_input <- paste(
  rep("-", times = banner_input_length - nchar("-- Input data ")),
  collapse = ''
 )


}
