library(tidyverse)
library(riskRegression)
library(survival)

tictoc::tic()
fit <- orsf(pbc_orsf, formula = time+status ~ . - id)
tictoc::toc()

sink("orsf-output.txt")
pd_vals <- orsf_pd_oob(fit,
                       expand_grid = FALSE,
                       pred_spec = list(bili = 1:5, trt = 'placebo'),
                       pred_horizon = seq(100, 1000, by=100))
sink()

fit$importance->tmp

tictoc::tic()
orsf_pd_oob(fit, pred_spec = list(bili = c(1:5)))
tictoc::toc()

prd_5 = predict(fit, new_data = pbc_orsf, n_thread = 5, pred_type = 'mort',
                pred_aggregate = F, pred_horizon = c(500, 1000))

microbenchmark::microbenchmark(
 prd_1 = predict(fit, new_data = pbc_orsf, n_thread = 1, pred_type = 'leaf',
                  pred_aggregate = F),

 prd_5 = predict(fit, new_data = pbc_orsf, n_thread = 5, pred_type = 'leaf',
                  pred_aggregate = F)
)

max(prd_1 - prd_5)


res <- oob <- vector(mode = 'numeric', length = 100)

for(i in seq(100)){

 train <- sample(nrow(pbc_orsf), 150)

 fit <- orsf(pbc_orsf[train, ],
             formula = Surv(time, status) ~ . - id,
             oobag_pred_type = 'surv',
             oobag_pred_horizon = 1000,
             split_rule = 'logrank',
             n_thread = 10)

 oob[i] = as.numeric(fit$eval_oobag$stat_values)

 prd <- predict(fit,
                new_data = pbc_orsf[-train, ],
                pred_horizon = 1000,
                pred_type = 'risk',
                n_thread = 10)

 y_mat <- as.matrix(pbc_orsf[-train, c('time', 'status')])
 w_vec <- rep(1, nrow(y_mat))
 s_vec <- 1-prd

 res[i] = oobag_c_survival(y_mat, w_vec, s_vec)


}

mean(oob)
mean(res)
mean(oob-res)

# sink()

library(randomForestSRC)

microbenchmark::microbenchmark(
 aorsf_5 = orsf(pbc_orsf, Surv(time, status) ~ . - id,
                n_tree = 500,
                mtry = 3,
                leaf_min_obs = 10,
                n_split = 5,
                importance = 'none',
                n_thread = 5),
 rfsrc = randomForestSRC::rfsrc(Surv(time, status) ~ . -id,
                                ntree = 500,
                                mtry = 3,
                                nthread = 5,
                                samptype = 'swr',
                                importance = 'none',
                                nodesize = 10,
                                nsplit = 5,
                                data = as.data.frame(pbc_orsf))
)

# sink()

fit$eval_oobag


fit$forest[-1] |>
 as_tibble() |>
 slice(1) |>
 unnest(everything()) |>
 mutate(node_id = seq(0, n()-1),
        .before = 1) |>
 print(n=100)

flchain_orsf <- flchain |>
 rename(time = futime, status = death) |>
 select(-chapter) |>
 filter(time > 0) %>%
 tidyr::drop_na()

microbenchmark::microbenchmark(

 aorsf_5 = orsf(flchain_orsf, Surv(time, status) ~ .,
                n_tree = 500,
                mtry = 3,
                leaf_min_obs = 10,
                n_split = 5,
                importance = 'permute',
                n_thread = 5),

 rfsrc = randomForestSRC::rfsrc(Surv(time, status) ~ .,
                                ntree = 500,
                                mtry = 3,
                                nthread = 5,
                                samptype = 'swr',
                                importance = 'permute',
                                nodesize = 10,
                                nsplit = 5,
                                data = as.data.frame(flchain_orsf)),

 times = 3

)


