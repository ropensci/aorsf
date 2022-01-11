
#' ORSF summary of univariate information
#'
#' ORSF's linear combinations of inputs can be used to provide helpful data
#'   about individual variables. Summarizing the univariate information from
#'   an ORSF provides data on the importance of individual variables and the
#'   expected predicted risk at designated values of the variables.
#'
#' @inheritParams predict.aorsf
#'
#' @param n_variables (_integer_) how many variables should be summarized?
#'   Setting this input to a lower number will improve computation time.
#'
#' @return an object of class 'aorsf_summary'
#'
#' @export
#'
#' @examples
#'
#' object <- orsf(pbc_orsf, Surv(time, status) ~ . - id)
#'
#' orsf_summarize_uni(object, n_variables = 3)
#'
orsf_summarize_uni <- function(object,
                               n_variables = NULL,
                               times = NULL,
                               risk = TRUE){

 # for CRAN check:
 medn <- name <- value <- level <- variable <- NULL

 check_arg_is(arg_value = object,
              arg_name = 'object',
              expected_class = 'aorsf')

 if(!is.null(n_variables)){

  check_arg_type(arg_value = n_variables,
                 arg_name = 'n_variables',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = n_variables,
                       arg_name = 'n_variables')

  check_arg_gteq(arg_value = n_variables,
                 arg_name = 'n_variables',
                 bound = 1)

  check_arg_lteq(arg_value = n_variables,
                 arg_name = 'n_variables',
                 bound = length(get_names_x(object, one_hot_names = FALSE)),
                 append_to_msg = "(total number of predictors)")


  check_arg_length(arg_value = n_variables,
                   arg_name = 'n_variables',
                   expected_length = 1)

 }

 if(!is.null(times)){

  check_arg_type(arg_value = times,
                 arg_name = 'times',
                 expected_type = 'numeric')

  check_arg_gt(arg_value = times,
               arg_name = 'times',
               bound = 0)

 }

 if(is.null(times)) times <- object$time_pred

 x_numeric_key <- get_numeric_bounds(object)

 fctr_info <- get_fctr_info(object)

 n_obs <- get_n_obs(object)

 importance <- orsf_vi(object, group_factors = TRUE)

 if(is.null(n_variables)) n_variables <- length(importance)

 pd_spec <- list_init(names(importance)[seq(n_variables)])

 for(x_name in names(pd_spec)){

  if(x_name %in% colnames(x_numeric_key)){

   pd_spec[[x_name]] <- unique(
    as.numeric(x_numeric_key[c('25%','50%','75%'), x_name])
   )

  } else if (x_name %in% fctr_info$cols) {

   pd_spec[[x_name]] <- fctr_info$lvls[[x_name]]

  }

 }

 pd_output <- orsf_pd_summary(object = object,
                              pd_spec = pd_spec,
                              expand_grid = FALSE,
                              risk = risk,
                              prob_values = c(0.25, 0.50, 0.75),
                              times = times)

 fctrs_unordered <- c()

 # did the orsf have factor variables?
 if(!is_empty(fctr_info$cols)){
  fctrs_unordered <- fctr_info$cols[!fctr_info$ordr]
 }

 # some cart-wheels here for backward compatibility.
 name_rep <- rle(as.integer(as.factor(pd_output$variable)))

 pd_output$importance <- rep(importance[name_rep$values],
                             times = name_rep$lengths)

 pd_output[, value := fifelse(test = is.na(value),
                              yes = level,
                              no = table.glue::table_value(value))]

 # if a := is used inside a function with no DT[] before the end of the
 # function, then the next time DT or print(DT) is typed at the prompt,
 # nothing will be printed. A repeated DT or print(DT) will print.
 # To avoid this: include a DT[] after the last := in your function.
 pd_output[]

 setcolorder(pd_output, c('variable',
                          'importance',
                          'value',
                          'mean',
                          'medn',
                          'lwr',
                          'upr'))

 structure(
  .Data = list(dt = pd_output,
               risk = risk,
               times = times),
  class = 'aorsf_summary_uni'
 )


}

#' Print ORSF summary
#'
#' @param x an object of class 'aorsf_summary'
#'
#' @param n_variables The number of variables to print
#'
#' @param ... not used
#'
#' @return nothing - output is printed to console.
#'
#' @export
#'
#' @examples
#'
#' object <- orsf(pbc_orsf, Surv(time, status) ~ . - id)
#'
#' smry <- orsf_summarize_uni(object, n_variables = 3)
#'
#' print(smry)
#'
print.aorsf_summary_uni <- function(x, n_variables = NULL, ...){


 if(is.null(n_variables)) n_variables <- length(unique(x$dt$variable))

 risk_or_surv <- if(x$risk) "risk" else "survival"

 msg_btm <- paste("Predicted", risk_or_surv,
                  "at time t =", x$times,
                  "for top", n_variables,
                  "predictors")

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

 name_index <- rle(as.integer(as.factor(x$dt$variable)))

 row_current <- 1

 i_vals <- seq(min(n_variables, length(name_index$values)))

 for(i in i_vals){

  name <- paste0('-- ', name_index$values[i],
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
   (banner_input_length - banner_value_length - nchar(risk_or_surv)) / 2

  header_length <- header_length - 1.5

  header_row <- paste(
   paste(rep(" ", times = banner_value_length), collapse = ''),
   paste(c("|",rep("-", times = header_length)), collapse = ''),
   " ",
   risk_or_surv,
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

}
