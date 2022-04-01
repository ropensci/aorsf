
# NOTE: THIS IS STILL EXPERIMENTAL; I MAY REMOVE IT.
# (why? Using vip's interaction function works better, but is slower)

#' ORSF interactions
#'
#' @srrstats {G1.4} *documented with Roxygen*
#'
#' @inheritParams predict.aorsf
#'
#' @param min_pairwise_obs (_integer_) minimum number of observations
#'   where both variables were included in a linear combination together.
#'   The default is the number of trees in `object` divided by the mean
#'   number of leaves in the trees. Default is the number of trees
#'   divided by the average number of leaves per tree, rounded to the
#'   nearest integer, i.e.,
#'   `min_pairwise_obs = round(number of trees / mean(number of leaves))`.
#'
#' @return a `data.frame` with pairwise interaction scores for each
#'   pair of predictor variables in `object`.
#'
#'
#' @srrstats {G2.0a} *specifying the length of `min_pairwise_obs`.*
#'
#' @details `min_pairwise_obs` should be a single value and should be large
#'  enough to prevent consideration of variable pairs that are so infrequently
#'  used together that it would not make sense to consider them interacting.
#'  For example, if two variables are only used together in one split of
#'  an entire random forest (i.e., number of pairwise observations = 1),
#'  it would not make sense to compute a correlation.
#'
#' @examples
#'
#' set.seed(329)
#'
#' fit <- orsf(pbc_orsf, Surv(time, status) ~ . - id, n_tree = 2500)
#'
#' intr <- aorsf:::orsf_interaction(fit)
#'
#' # edema==1 and bili are strongest interacting pair
#' print(intr)
#'
#' # make a list containing the variable values you
#' # want to compute partial dependence for
#' pd_spec <- list(edema = c("0", "0.5", "1"),
#'                 bili = seq(0.6, 7.1, by = 0.5))
#'
#' # orsf_pd_summary automatically computes pd for all combinations
#' # in the list (this can be turned off with expand_grid = FALSE)
#' pd_data <- orsf_pd_summary(object = fit,
#'                            pd_spec = pd_spec,
#'                            expand_grid = TRUE)
#'
#' # aligning predictions at lowest value of bili
#' min_ed_0 <- with(pd_data, mean[edema == "0"   & bili == 0.6])
#' min_ed_1 <- with(pd_data, mean[edema == "0.5" & bili == 0.6])
#' min_ed_2 <- with(pd_data, mean[edema == "1"   & bili == 0.6])
#'
#' pd_data_aligned <-
#'  within(pd_data, {
#'   value <- mean
#'   value[edema == "0" ] <- value[edema == "0" ] - min_ed_0
#'   value[edema == "0.5"] <- value[edema == "0.5"] - min_ed_1
#'   value[edema == "1" ] <- value[edema == "1" ] - min_ed_2
#'  })
#'
#' head(pd_data_aligned)
#'
#' library(ggplot2)
#'
#' ggplot(pd_data_aligned) +
#'  aes(x = bili, y = value, col = edema) +
#'  geom_line() +
#'  labs(y = 'Predictions centered at Bilirubin = 0.6',
#'       x = 'Bilirubin',
#'       title = 'Interaction between bilirubin and edema')


orsf_interaction <- function(object,
                             min_pairwise_obs = NULL){

 #' @srrstats {G2.8} *As part of initial pre-processing, run checks on inputs to ensure that all other sub-functions receive inputs of a single defined class or type.*

 check_arg_is(arg_value = object,
              arg_name = "object",
              expected_class = 'aorsf')

 if(!is.null(min_pairwise_obs)){

  check_arg_type(arg_value = min_pairwise_obs,
                 arg_name = 'min_pairwise_obs',
                 expected_type = 'numeric')

  check_arg_is_integer(arg_value = min_pairwise_obs,
                       arg_name = 'min_pairwise_obs')

  check_arg_gteq(arg_value = min_pairwise_obs,
                 arg_name = 'min_pairwise_obs',
                 bound = 1)

  check_arg_length(arg_value = min_pairwise_obs,
                   arg_name = 'min_pairwise_obs',
                   expected_length = 1)

 }

 if(is.null(min_pairwise_obs))
  min_pairwise_obs <- round(get_n_tree(object) / get_n_leaves_mean(object))

 # for CRAN:
 cor <- value <- NULL

 dt_betas <- dt_means <- list()

 xnames <- get_names_x(object, ref_code_names = TRUE)

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

  # imat[xnames[-i], xnames[i]] <-
  #  t(cor(dt_betas[,i], dt_means[,-i], use = 'pairwise.complete.obs')^2)
  # tmp = vector(mode = 'numeric', length = n_vars)

  observed_i <- !is.na(dt_betas[,i])

  for(j in setdiff(seq(n_vars), i) ){

   observed_j <- !is.na(dt_means[,j])

   pairwise_complete <- observed_i & observed_j

   n_pairwise_complete <- sum(pairwise_complete)

   if(n_pairwise_complete > min_pairwise_obs){

    imat[i, j] <-
     cor(dt_betas[pairwise_complete, i],
         dt_means[pairwise_complete, j])^2
   }

  }


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
