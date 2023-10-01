library(tidyverse)
library(riskRegression)
library(survival)

tictoc::tic()
fit <- orsf(pbc_orsf,
            formula = time+status ~ . - id,
            oobag_pred_type = 'none',
            sample_with_replacement = FALSE,
            sample_fraction = 2/3)
tictoc::toc()

all(fit$data == pbc_orsf)

tmp <- as.data.frame(cbind(y=as.numeric(fit$pred_oobag), x=pbc_orsf$bili))

plot(x=tmp$x, y=tmp$y)

sink("orsf-output.txt")
pd_vals <- orsf_pd_oob(fit,
                       expand_grid = FALSE,
                       pred_type = 'risk',
                       pred_spec = list(bili = 1:5,
                                        sex = c("m", "f")),
                       pred_horizon = c(1000, 2000, 4000))
sink()

fit$importance->tmp

tictoc::tic()
orsf_pd_oob(fit, pred_spec = list(bili = c(1:5)))
tictoc::toc()

rfsrc_fit = randomForestSRC::rfsrc(Surv(time, status) ~ . -id,
                               ntree = 500,
                               mtry = 3,
                               nthread = 5,
                               samptype = 'swr',
                               importance = 'none',
                               nodesize = 10,
                               nsplit = 5,
                               data = as.data.frame(pbc_orsf))

microbenchmark::microbenchmark(
 prd_aorsf = predict(fit, new_data = pbc_orsf, n_thread = 10),
 prd_rfsrc = predict(rfsrc_fit, newdata = pbc_orsf)
)

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
                                data = as.data.frame(pbc_orsf)),
 times = 5
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
                importance = 'none',
                n_thread = 5),

 rfsrc = randomForestSRC::rfsrc(Surv(time, status) ~ .,
                                ntree = 500,
                                mtry = 3,
                                nthread = 5,
                                samptype = 'swr',
                                importance = 'none',
                                nodesize = 10,
                                nsplit = 5,
                                data = as.data.frame(flchain_orsf)),

 times = 3

)


