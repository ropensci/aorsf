
#' General rounding for presentation
#'
#'  casts numeric vectors into character vectors for presentation.
#'
#' @param x a vector of numeric values.
#'
#' @return a vector of character values (rounded numbers).
#'
#' @noRd
#'

round_magnitude <- function(x){

 out <- rep(NA_character_, length(x))

 if (all(is.na(x))) return(out)

 # take absolute value to round based on magnitude
 x_abs <- abs(x)

 breaks <- c(0, 1, 10, Inf)
 decimals <- c(2, 1, 0)

 # x_cuts create boundary categories for rounding
 x_cuts <- cut(
  x_abs,
  breaks = breaks,
  include.lowest = TRUE,
  right = FALSE
 )

 out_breaks <- lapply(
  levels(x_cuts),
  function(.x) which(x_cuts == .x)
 )

 for (i in seq_along(out_breaks)) {

  ob <- out_breaks[[i]]

  if(!is_empty(ob)) {

   ob_rounded <- round(x[ob], digits = decimals[i])

   out[ob] <- format(ob_rounded, nsmall = decimals[i], trim = TRUE)

  }

 }

 out

}
