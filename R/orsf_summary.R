
#' Summarize ORSF model
#'
#' @inheritParams predict.aorsf
#' @param n_variables how many variables should be summarized?
#' @return an object of class 'aorsf_summary'
#' @export
#'
#' @examples
#'
#' object <- orsf(pbc_orsf, Surv(time, status) ~ . - id)
#'
#' orsf_summarize_uni(object)
#'
orsf_summarize_uni <- function(object,
                               n_variables = NULL,
                               times = NULL,
                               risk = TRUE){

 # for CRAN check:
 medn <- name <- value <- level <- variable <- NULL

 if(is.null(times))
  times <- object$time_pred

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


 # TODO: write this in cpp?

 # position <- 1
 #
 # pratio <- pdiff <- rep(NA_real_, nrow(pd_output))
 #
 # pred <- pd_output$pred
 # name <- pd_output$name
 #
 # while(position < nrow(pd_output)){
 #
 #  ref_index <- seq(position, position + n_obs - 1)
 #
 #  step <- 1
 #  pratio[ref_index] <- 1
 #  pdiff[ref_index] <- 0
 #
 #  new_index <- ref_index + n_obs
 #
 #  while( name[ new_index[1] ] == name[ ref_index[1] ] ){
 #
 #   step <- step + 1
 #   pratio[new_index] <- pred[new_index] / pred[ref_index]
 #   pdiff[new_index] <- pred[new_index] - pred[ref_index]
 #   new_index <- ref_index + n_obs * step
 #
 #   if(new_index[1] > nrow(pd_output)) break
 #
 #  }
 #
 #  position <- new_index[1]
 #
 # }
 #
 # pd_output$pred_ratio <- pratio
 # pd_output$pred_diff <- pdiff

 # end TODO

 fctrs_unordered <- c()

 if(!is_empty(fctr_info$cols)){
  fctrs_unordered <- fctr_info$cols[!fctr_info$ordr]
 }

 name_rep <- rle(pd_output$name)

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

 setcolorder(pd_output, c('name',
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

as.data.table.aorsf_summary_uni <- function(x) x$dt

#' Print ORSF summary
#'
#' @param x an object of class 'aorsf_summary'
#' @param n_variables The number of variables to print
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
#' summary(object)
#'
print.aorsf_summary_uni <- function(x, n_variables = 3, ...){


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

 .sd_fncy <- c("Variable Value",
               paste(c("Mean", "Median"), risk_or_surv),
               "25th Percentile",
               "75th Percentile")

 setnames(x$dt,
          old = .sd_orig,
          new = .sd_fncy)

 banner_input_length <-
  vapply(
   utils::capture.output(
    do.call(
     print,
     list(x = x$dt[1L, .SD, .SDcols = .sd_fncy],
          trunc.cols = TRUE)
    )
   ),
   nchar,
   integer(1)
  )[1L]

 name_index <- rle(x$dt$name)
 row_current <- 1

 i_vals <- seq(min(n_variables, length(name_index$values)))

 for(i in i_vals){

  name <- paste0('-- ', name_index$values[i],
                 " (VI Rank: ", i, ")")

  banner_input <- paste(
   rep("-", times = max(0, banner_input_length - nchar(name))),
   collapse = ''
  )

  cat("\n",
      name,
      " ",
      banner_input,
      "\n\n",
      sep = "")

  row_new <- row_current + name_index$lengths[i]-1

  print(x$dt[row_new:row_current, .SD, .SDcols = .sd_fncy],
        row.names = F,
        col.names = "top",
        trunc.cols = TRUE,
        digits = 4)

  row_current <- row_new+1

 }

 cat("\n", msg_btm)

 setnames(x$dt,
          old = .sd_fncy,
          new = .sd_orig)

}
