

#' ORSF interactions
#'
#' @inheritParams predict.aorsf
#'
#' @return a `data.frame` with pairwise interaction scores for each
#'   pair of predictor variables in `object`.
#'
#' @export
#'
#' @examples
#'
#' pbc_orsf$ascites <- factor(pbc_orsf$ascites)
#'
#' set.seed(32987)
#'
#' fit <- orsf(pbc_orsf, Surv(time, status) ~ . - id)
#'
#' intr <- orsf_interaction(fit)
#'
#' pd_spec <- list(ascites = c("0","1"),
#'                 bili = seq(0.6, 7.1, by = 0.25))
#'
#' pd_data <- orsf_pd_summary(fit, pd_spec = pd_spec, times = 1000)
#'
#' # aligning predictions at lowest value of bili
#' min_asc_0 <- with(pd_data, mean[ascites == 0 & bili == 0.6])
#' min_asc_1 <- with(pd_data, mean[ascites == 1 & bili == 0.6])
#'
#' pd_data_aligned <-
#'  within(pd_data, {
#'   mean[ascites == 0] <- mean[ascites == 0] - min_asc_0
#'   mean[ascites == 1] <- mean[ascites == 1] - min_asc_1
#'  })
#'
#' library(ggplot2)
#'
#' ggplot(pd_data_aligned) +
#'  aes(x = bili, y = mean, col = factor(ascites)) +
#'  geom_line()

orsf_interaction <- function(object){

 check_arg_is(arg_value = object,
              arg_name = "object",
              expected_class = 'aorsf')

 # for CRAN:
 cor <- value <- NULL

 dt_betas <- dt_means <- list()

 xnames <- get_names_x(object, one_hot_names = TRUE)

 for(t in seq_along(object$forest)){

  tree <- object$forest[[t]]

  nodes <- which(tree$children_left > 0)

  tree_means <- tree_coefs <- matrix(NA_real_,
                                     nrow = length(nodes),
                                     ncol = length(xnames))

  colnames(tree_coefs) <- colnames(tree_means) <- xnames

  for(i in seq(nrow(tree_coefs))){

   ii <- nodes[i]

   cols <- tree$col_indices[, ii] + 1

   beta_nonzero <- which(tree$betas[, ii] != 0)

   tree_coefs[i, cols[beta_nonzero]] <- tree$betas[beta_nonzero, ii]
   tree_means[i, cols[beta_nonzero]] <- tree$x_mean[beta_nonzero, ii]

  }

  dt_betas[[t]] <- as.data.table(tree_coefs)
  dt_means[[t]] <- as.data.table(tree_means)

 }

 dt_betas <- as.matrix(rbindlist(dt_betas))
 dt_means <- as.matrix(rbindlist(dt_means))

 imat <- matrix(NA_real_,
                nrow = length(xnames),
                ncol = length(xnames))

 colnames(imat) <- rownames(imat) <- xnames

 n_vars <- length(xnames)

 i_vals <- seq(n_vars)

 for(i in i_vals){

  imat[xnames[-i], xnames[i]] <-
   t(cor(dt_betas[,i], dt_means[,-i], use = 'pairwise.complete.obs')^2)

 }

 for(i in i_vals){

  for(j in seq(i, n_vars)){
   if(j > i){
    imat[i, j] <- (imat[i, j] + imat[j, i]) / 2
    imat[j, i] <- NA_real_
   }

  }

 }

 dt <- as.data.table(imat, keep.rownames = 'v1')
 dt <- melt(dt, id.vars = 'v1', variable.name = 'v2', na.rm = TRUE)
 dt <- dt[order(value, decreasing = TRUE), ]
 # dt[]

 dt

}
