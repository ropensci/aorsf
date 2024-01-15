
#' Univariate summary
#'
#' @description Summarize the univariate information from an ORSF object
#'
#' @inheritParams orsf_pd_oob
#'
#' @param n_variables (_integer_) how many variables should be summarized?
#'   Setting this input to a lower number will reduce computation time.
#'
#' @param importance `r roxy_importance_header()`
#' - `r roxy_importance_none()`
#' - `r roxy_importance_anova()`
#' - `r roxy_importance_negate()`
#' - `r roxy_importance_permute()`
#'
#' For details on these methods, see [orsf_vi].
#'
#' @return an object of class 'orsf_summary', which includes data on
#'
#' - importance of individual predictors.
#' - expected values of predictions at specific values of predictors.
#'
#' @details
#'
#'  If `pred_horizon` is left unspecified, the median value of
#'    the time-to-event variable in `object`'s training data will be used.
#'    It is recommended to always specify your own prediction horizon,
#'    as the median time may not be an especially meaningful horizon to
#'    compute predicted risk values at.
#'
#'  If `object` already has variable importance values, you can
#'    safely bypass the computation of variable importance in this function
#'    by setting importance = 'none'.
#'
#' @export
#'
#' @seealso as.data.table.orsf_summary_uni
#'
#' @examples
#'
#' object <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 25)
#'
#' # since anova importance was used to make object, it is also
#' # used for ranking variables in the summary, unless we specify
#' # a different type of importance
#'
#' orsf_summarize_uni(object, n_variables = 2)
#'
#' # if we want to summarize object according to variables
#' # ranked by negation importance, we can compute negation
#' # importance within orsf_summarize_uni() as follows:
#'
#' orsf_summarize_uni(object, n_variables = 2, importance = 'negate')
#'
#'
orsf_summarize_uni <- function(object,
                               n_variables = NULL,
                               pred_horizon = NULL,
                               pred_type = NULL,
                               importance = NULL,
                               verbose_progress = FALSE,
                               ...){

 check_dots(list(...), .f = orsf_summarize_uni)

 check_arg_is(arg_value = object,
              arg_name = 'object',
              expected_class = 'ObliqueForest')

 object$summarize_uni(n_variables     = n_variables,
                      pred_horizon    = pred_horizon,
                      pred_type       = pred_type,
                      importance_type = importance,
                      verbose_progress = verbose_progress)

}

#' Print ORSF summary
#'
#' @param x an object of class 'orsf_summary'
#'
#' @param n_variables The number of variables to print
#'
#' @param ... `r roxy_dots()`
#'
#' @return invisibly, `x`
#'
#' @export
#'
#' @examples
#'
#' object <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 25)
#'
#' smry <- orsf_summarize_uni(object, n_variables = 2)
#'
#' print(smry)
#'
print.orsf_summary_uni <- function(x, n_variables = NULL, ...){

 check_dots(list(...), .f = print.orsf_summary_uni)

 if(is.null(n_variables)) n_variables <- length(unique(x$dt$variable))

 pred_label <- switch(
  x$pred_type,
  'risk' = 'Risk',
  'surv' = 'Survival probability',
  'chf'  = 'Cumulative hazard',
  'mort' = 'Mortality',
  'prob' = "Probability",
  "time" = "Survival time",
  "mean" = "Expected value"
 )

 extra_surv_text <- ""

 if(!is.null(x$pred_horizon))
    extra_surv_text <- paste0("at time t = ", x$pred_horizon, " ")

 msg_btm <- paste0(
  "Predicted ", tolower(pred_label), " ", extra_surv_text,
  "for top ", n_variables, " predictors"
 )

 .sd_orig <- c("value",
               "mean",
               "medn",
               "lwr",
               "upr")

 .sd_fncy <- c("Value",
               "Mean",
               "Median",
               "25th %",
               "75th %")

 if(any(.sd_fncy %in% names(x$dt)))
  setnames(x$dt, old = .sd_fncy, new = .sd_orig, skip_absent = TRUE)

 setnames(x$dt,
          old = .sd_orig,
          new = .sd_fncy)

 banner_input_length <-
  as.integer(
   vapply(
    utils::capture.output(
     do.call(
      print,
      list(x = x$dt[, .SD, .SDcols = .sd_fncy],
           row.names = FALSE,
           trunc.cols = TRUE)
     )
    ),
    nchar,
    integer(1)
   )[1L]
  )

 first_col_vals <- nchar(x$dt$Value)
 first_col_width <- max(first_col_vals)
 first_col_pad <- first_col_width - first_col_vals

 for(i in seq_along(x$dt$Value)){

  pad <- paste(
   rep(" ", times = first_col_pad[i]),
   collapse = ''
  )

  x$dt$Value[i] <- paste0(pad, x$dt$Value[i])

 }

 # cart-wheels for rle backward compatibility
 f <- as.factor(x$dt$variable)

 name_index <- rle(as.integer(f))

 row_current <- 1

 i_vals <- seq(min(n_variables, length(name_index$values)))

 for(i in i_vals){

  name <- paste0('-- ', levels(f)[name_index$values[i]],
                 " (VI Rank: ", i, ")")

  banner_input <- paste(
   rep("-", times = max(0, banner_input_length - nchar(name))),
   collapse = ''
  )

  banner_value_length <-
   as.integer(
    vapply(
     utils::capture.output(
      do.call(
       print,
       list(x = x$dt[i, .SD, .SDcols = .sd_fncy[1]],
            row.names = FALSE,
            trunc.cols = TRUE)
      )
     ),
     nchar,
     integer(1)
    )[1L]
   )

  banner_value_length <- banner_value_length + 1

  header_length <-
   (banner_input_length - banner_value_length - nchar(pred_label)) / 2

  header_length <- header_length - 1.5

  header_row <- paste(
   paste(rep(" ", times = banner_value_length), collapse = ''),
   paste(c("|",rep("-", times = header_length)), collapse = ''),
   " ",
   pred_label,
   " ",
   paste(c(rep("-", times = header_length), "|"), collapse = ''),
   collapse = '',
   sep = ''
  )

  cat("\n",
      name,
      " ",
      banner_input,
      "\n\n",
      header_row,
      "\n",
      sep = "")

  row_new <- row_current + name_index$lengths[i]-1

  print(x$dt[row_current:row_new, .SD, .SDcols = .sd_fncy],
        row.names = FALSE,
        col.names = "top",
        trunc.cols = TRUE)

  row_current <- row_new+1

 }

 cat("\n", msg_btm, "\n")

 setnames(x$dt,
          old = .sd_fncy,
          new = .sd_orig)

 invisible(x)

}



#' Coerce to data.table
#'
#' Convert an 'orsf_summary' object into a `data.table` object.
#'
#' @param x an object of class 'orsf_summary_uni'
#'
#' @param ... not used
#'
#' @return a [data.table][data.table::data.table-package]
#'
#' @export
#'
#' @examples
#'
#'\dontrun{
#'
#' library(data.table)
#'
#' object <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 25)
#'
#' smry <- orsf_summarize_uni(object, n_variables = 2)
#'
#' as.data.table(smry)
#'
#'}
#'
#'
as.data.table.orsf_summary_uni <- function(x, ...){
 x$dt
}





