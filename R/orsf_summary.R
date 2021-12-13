

#' Summarize ORSF model
#'
#' @inheritParams predict.aorsf
#'
#' @return an object of class 'aorsf_summary'
#' @export
#'
#' @examples
#'
#' object <- orsf(pbc_orsf, Surv(time, status) ~ . - id)
#'
#' summary(object)
#'
summary.aorsf <- function(object, times = NULL, risk = TRUE, ...){

 if(is.null(times))
  times <- object$time_pred

 x_numeric_key <- get_numeric_bounds(object)

 fctr_info <- get_fctr_info(object)

 n_obs <- get_n_obs(object)

 pd_spec <- list_init(get_names_x(object))

 for(x_name in names(pd_spec)){

  if(x_name %in% colnames(x_numeric_key)){

   pd_spec[[x_name]] <- unique(
    as.numeric(x_numeric_key[c('25%','50%','75%'), x_name])
   )

  } else if (x_name %in% fctr_info$cols) {

   pd_spec[[x_name]] <- fctr_info$lvls[[x_name]]

  }

 }

 pd_data <- orsf_pd_ice(object = object,
                        pd_spec = pd_spec,
                        expand_grid = FALSE,
                        risk = risk,
                        times = times)



 importance <- orsf_vi(object)

 # TODO: write this in cpp?

 # position <- 1
 #
 # pratio <- pdiff <- rep(NA_real_, nrow(pd_data))
 #
 # pred <- pd_data$pred
 # name <- pd_data$name
 #
 # while(position < nrow(pd_data)){
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
 #   if(new_index[1] > nrow(pd_data)) break
 #
 #  }
 #
 #  position <- new_index[1]
 #
 # }
 #
 # pd_data$pred_ratio <- pratio
 # pd_data$pred_diff <- pdiff

 # end TODO

 pd_smry <- pd_data[, .(pred_mean   = mean(pred),
                        pred_median = median(pred)),
                    by = list(name, value, level)]

 fctrs_unordered <- c()

 if(!is_empty(fctr_info$cols)){
  fctrs_unordered <- fctr_info$cols[!fctr_info$ordr]
 }


 pd_smry[, variable := fifelse(test = name %in% fctrs_unordered,
                               yes = paste(name, level, sep = '_'),
                               no = name)]


 # lvls_percentile <- paste(c('25th','50th','75th'), 'percentile')
 #
 # pd_smry[, level := fifelse(test = is.na(level),
 #                            yes = lvls_percentile[seq(.N)],
 #                            no = level),
 #         by = name]


 vi_smry <- data.table(variable = names(importance),
                       importance = as.numeric(importance))

 type_smry <- data.table(name = get_names_x(object),
                         type = get_types_x(object))

 out <- vi_smry[pd_smry, on = 'variable']

 out[, importance := mean(importance, na.rm=TRUE), by = name]

 out <- out[type_smry, on = 'name']

 setorder(out, -importance)

 out[, variable := NULL]

 # if a := is used inside a function with no DT[] before the end of the
 # function, then the next time DT or print(DT) is typed at the prompt,
 # nothing will be printed. A repeated DT or print(DT) will print.
 # To avoid this: include a DT[] after the last := in your function.
 out[]


 setcolorder(out, c('name',
                    'importance',
                    'type',
                    'value',
                    'level',
                    'pred_mean',
                    'pred_median'))
 out


}
