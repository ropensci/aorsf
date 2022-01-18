



#' ORSF presentation
#'
#' @srrstats {G1.4} *documented with Roxygen*
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

 info_oobag_c <- ifelse(test = contains_oobag(x),
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
