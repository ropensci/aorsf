

get_max_time        <- function(object) attr(object, 'max_time')
get_n_leaves_mean   <- function(object) attr(object, 'n_leaves_mean')
get_mtry            <- function(object) attr(object, 'mtry')
get_max_time        <- function(object) attr(object, 'max_time')
get_n_obs           <- function(object) attr(object, 'n_obs')
get_n_tree          <- function(object) attr(object, 'n_tree')
get_n_events        <- function(object) attr(object, 'n_events')
get_fctr_info       <- function(object) attr(object, 'fctr_info')
get_names_y         <- function(object) attr(object, 'names_y')
get_names_x         <- function(object) attr(object, 'names_x')
get_types_x         <- function(object) attr(object, 'types_x')
get_leaf_min_events <- function(object) attr(object, 'leaf_min_events')
get_leaf_min_obs    <- function(object) attr(object, 'leaf_min_obs')
get_numeric_bounds  <- function(object) attr(object, 'numeric_bounds')

has_oobag <- function(object) !is_empty(object$surv_oobag)
