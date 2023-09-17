library(survival)
library(tidyverse)
library(riskRegression)


.pbc_orsf <- pbc_orsf %>%
 mutate(stage = factor(stage, ordered = F))

x <- model.matrix(~. -1, data = select(.pbc_orsf, -time, -status, -id))
y <- as.matrix(.pbc_orsf[, c('time', 'status')])

# .flchain <- flchain |>
#  rename(time = futime, status = death) |>
#  select(-chapter) |>
#  tidyr::drop_na()
# x <- model.matrix(~. -1, data = select(.flchain, -time, -status))
# y <- as.matrix(.flchain[, c('time', 'status')])

w <- rep(1, nrow(x))

sorted <-
 collapse::radixorder(y[, 1],  # order this way for risk sets
                      -y[, 2]) # order this way for oob C-statistic.

y <- y[sorted, ]
x <- x[sorted, ]
w <- w[sorted]

f <- function(x, y, w){
 matrix(runif(ncol(x)), ncol=1)
}

pred_horizon <- median(y[, 'time'])

sink("orsf-output.txt")

orsf_tree = aorsf:::orsf_cpp(x,
                             y,
                             w,
                             tree_type_R = 3,
                             tree_seeds = 1:500,
                             loaded_forest = list(),
                             n_tree = 500,
                             mtry = 3,
                             vi_type_R = 2,
                             vi_max_pvalue = 0.01,
                             lincomb_R_function = f,
                             oobag_R_function = f,
                             leaf_min_events = 5,
                             leaf_min_obs = 5,
                             split_rule_R = 1,
                             split_min_events = 5,
                             split_min_obs = 10,
                             split_min_stat = 0,
                             split_max_cuts = 5,
                             split_max_retry = 3,
                             lincomb_type_R = 1,
                             lincomb_eps = 1e-9,
                             lincomb_iter_max = 1,
                             lincomb_scale = TRUE,
                             lincomb_alpha = 1,
                             lincomb_df_target = 1,
                             lincomb_ties_method = 1,
                             pred_type_R = 1,
                             pred_mode = FALSE,
                             pred_horizon = pred_horizon,
                             oobag = TRUE,
                             oobag_eval_every = 500,
                             n_thread = 1)

sink()

orsf_tree$forest[-1] |>
 as_tibble() |>
 slice(1) |>
 unnest(everything()) |>
 mutate(node_id = seq(0, n()-1),
        .before = 1) |>
 print(n=100)


# sink("orsf-output.txt")
res_cstat_oobag <- res_ipa_oobag <- res_cstat_new <- res_ipa_new <- c()

n_train <- nrow(pbc_orsf)/2

for(i in seq(5)){

 print(i)

 train_rows <- sort(sample(nrow(x), size = n_train))

 orsf_fit_1 = aorsf:::orsf_cpp(x[train_rows, ],
                               y[train_rows, ],
                               w[train_rows],
                               tree_seeds = 1:500,
                               loaded_forest = list(),
                               n_tree = 500,
                               mtry = 5,
                               vi_type_R = 3,
                               vi_max_pvalue = 0.01,
                               lincomb_R_function = f,
                               oobag_R_function = f,
                               leaf_min_events = 1,
                               leaf_min_obs = 5,
                               split_rule_R = 1,
                               split_min_events = 5,
                               split_min_obs = 10,
                               split_min_stat = 0,
                               split_max_cuts = 5,
                               split_max_retry = 3,
                               lincomb_type_R = 1,
                               lincomb_eps = 1e-9,
                               lincomb_iter_max = 1,
                               lincomb_scale = TRUE,
                               lincomb_alpha = 1,
                               lincomb_df_target = 1,
                               lincomb_ties_method = 1,
                               pred_type_R = 1,
                               pred_mode = FALSE,
                               pred_horizon = pred_horizon,
                               oobag = T,
                               oobag_eval_every = 500,
                               n_thread = 1)

 orsf_fit_5 = aorsf:::orsf_cpp(x[train_rows, ],
                               y[train_rows, ],
                               w[train_rows],
                               tree_seeds = 1:500,
                               loaded_forest = list(),
                               n_tree = 500,
                               mtry = 5,
                               vi_type_R = 3,
                               vi_max_pvalue = 0.01,
                               lincomb_R_function = f,
                               oobag_R_function = f,
                               leaf_min_events = 1,
                               leaf_min_obs = 5,
                               split_rule_R = 1,
                               split_min_events = 5,
                               split_min_obs = 10,
                               split_min_stat = 0,
                               split_max_cuts = 5,
                               split_max_retry = 3,
                               lincomb_type_R = 1,
                               lincomb_eps = 1e-9,
                               lincomb_iter_max = 1,
                               lincomb_scale = TRUE,
                               lincomb_alpha = 1,
                               lincomb_df_target = 1,
                               lincomb_ties_method = 1,
                               pred_type_R = 1,
                               pred_mode = FALSE,
                               pred_horizon = pred_horizon,
                               oobag = T,
                               oobag_eval_every = 500,
                               n_thread = 10)


 expect_equal(orsf_fit_1$vi_numer, orsf_fit_5$vi_numer)
 expect_equal(orsf_fit_1$vi_denom, orsf_fit_5$vi_denom)

 expect_equal(orsf_fit_1$forest, orsf_fit_5$forest)
 expect_equal(orsf_fit_1$pred_oobag, orsf_fit_5$pred_oobag)

 # sink()
 #
 # data.frame(variable = colnames(x),
 #            vi = orsf_fit$vi_numer/orsf_fit$vi_denom) |>
 #  dplyr::arrange(vi)
 #
 orsf_fit_1$forest[-1] |>
  as_tibble() |>
  slice(1) |>
  unnest(everything()) |>
  mutate(node_id = seq(0, n()-1),
         .before = 1) |>
  print(n=100)
 #
 # sink("orsf-output.txt")

 orsf_pred_1 = aorsf:::orsf_cpp(x[-train_rows, ],
                                y[-train_rows, ],
                                w[-train_rows],
                                tree_seeds = 1:500,
                                loaded_forest = orsf_fit_1$forest,
                                n_tree = 500,
                                mtry = 5,
                                vi_type_R = 3,
                                vi_max_pvalue = 0.01,
                                lincomb_R_function = f,
                                oobag_R_function = f,
                                leaf_min_events = 1,
                                leaf_min_obs = 5,
                                split_rule_R = 1,
                                split_min_events = 10,
                                split_min_obs = 20,
                                split_min_stat = 3.84,
                                split_max_cuts = 5,
                                split_max_retry = 3,
                                lincomb_type_R = 1,
                                lincomb_eps = 1e-9,
                                lincomb_iter_max = 1,
                                lincomb_scale = TRUE,
                                lincomb_alpha = 1,
                                lincomb_df_target = 1,
                                lincomb_ties_method = 1,
                                pred_type_R = 1,
                                pred_mode = TRUE,
                                pred_horizon = c(pred_horizon),
                                oobag = FALSE,
                                oobag_eval_every = 500,
                                n_thread = 1)

 # sink("orsf-output.txt")

 orsf_pred_10 = aorsf:::orsf_cpp(x[-train_rows, ],
                                 y[-train_rows, ],
                                 w[-train_rows],
                                 tree_seeds = 1:500,
                                 loaded_forest = orsf_fit_1$forest,
                                 n_tree = 500,
                                 mtry = 5,
                                 vi_type_R = 3,
                                 vi_max_pvalue = 0.01,
                                 lincomb_R_function = f,
                                 oobag_R_function = f,
                                 leaf_min_events = 1,
                                 leaf_min_obs = 5,
                                 split_rule_R = 1,
                                 split_min_events = 10,
                                 split_min_obs = 20,
                                 split_min_stat = 3.84,
                                 split_max_cuts = 5,
                                 split_max_retry = 3,
                                 lincomb_type_R = 1,
                                 lincomb_eps = 1e-9,
                                 lincomb_iter_max = 1,
                                 lincomb_scale = TRUE,
                                 lincomb_alpha = 1,
                                 lincomb_df_target = 1,
                                 lincomb_ties_method = 1,
                                 pred_type_R = 1,
                                 pred_mode = TRUE,
                                 pred_horizon = c(pred_horizon),
                                 oobag = FALSE,
                                 oobag_eval_every = 500,
                                 n_thread = 5)

 # sink()

 pred_list <- list(one_thread = 1-orsf_pred_1$predictions,
                   ten_thread = 1-orsf_pred_10$predictions)

 expect_equal(pred_list$one_thread, pred_list$ten_thread)

 sc_oobag = Score(object = list(1-orsf_fit_1$pred_oobag),
                  Surv(time, status) ~ 1,
                  data = as.data.frame(y)[train_rows, ],
                  times = pred_horizon,
                  summary = 'IPA')

 sc = Score(object = pred_list,
            Surv(time, status) ~ 1,
            data = as.data.frame(y)[-train_rows, ],
            times = pred_horizon,
            summary = 'IPA')

 res_cstat_oobag <- c(res_cstat_oobag, sc_oobag$AUC$score$AUC[1])
 res_ipa_oobag <- c(res_ipa_oobag, sc_oobag$Brier$score$IPA[2])
 res_cstat_new <- c(res_cstat_new, sc$AUC$score$AUC[1])
 res_ipa_new <- c(res_ipa_new, sc$Brier$score$IPA[2])

}

mean(res_cstat_oobag)
mean(res_cstat_new)
mean(res_ipa_oobag)
mean(res_ipa_new)

microbenchmark::microbenchmark(

 orsf_fit = aorsf:::orsf_cpp(x,
                             y,
                             w,
                             tree_seeds = c(1:500),
                             loaded_forest = list(),
                             n_tree = 500,
                             mtry = 3,
                             vi_type_R = 2,
                             vi_max_pvalue = 0.01,
                             lincomb_R_function = f,
                             oobag_R_function = f,
                             leaf_min_events = 5,
                             leaf_min_obs = 10,
                             split_rule_R = 1,
                             split_min_events = 10,
                             split_min_obs = 20,
                             split_min_stat = 3.84,
                             split_max_cuts = 5,
                             split_max_retry = 3,
                             lincomb_type_R = 1,
                             lincomb_eps = 1e-5,
                             lincomb_iter_max = 1,
                             lincomb_scale = TRUE,
                             lincomb_alpha = 1,
                             lincomb_df_target = 1,
                             lincomb_ties_method = 1,
                             pred_type_R = 1,
                             pred_mode = FALSE,
                             pred_horizon = pred_horizon,
                             oobag = TRUE,
                             oobag_eval_every = 500,
                             n_thread = 10),

 # ranger_fit = ranger::ranger(x = x, y = y,
 #                             mtry = 3,
 #                             num.threads = 10,
 #                             num.trees = 500),

 rfsrc_fit = randomForestSRC::rfsrc(Surv(time, status) ~ .,
                                    ntree = 500,
                                    mtry = 3,
                                    nthread = 10,
                                    samptype = 'swr',
                                    importance = 'permute',
                                    nodesize = 10,
                                    data = as.data.frame(cbind(y, x))),

 times = 3


)

