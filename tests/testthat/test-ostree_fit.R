
library(survival)
library(randomForestSRC)
library(riskRegression)
devtools::load_all()

.pbc <-  pbc[order(pbc$time), ]
.pbc <- .pbc[complete.cases(.pbc), ]
.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1

x <- as.matrix(.pbc[, c('trt','age','ascites','hepato',
                        'spiders','chol','edema','bili',
                        'albumin','copper','alk.phos',
                        'ast','trig','platelet','protime',
                        'stage')])

y <- Surv(.pbc$time, .pbc$status)

set.seed(1)

microbenchmark::microbenchmark(
 tmp = orsf_fit(x,
                y,
                n_tree = 50,
                mtry_ = 3,
                leaf_min_events_ = 5,
                leaf_min_obs_ = 50,
                oobag_pred_ = TRUE,
                cph_iter_max_ = 10,
                cph_pval_max_ = .8)
)




res <- list()

for(i in 1:100){

 train <- sort(sample(nrow(x), 200))

 # set.seed(2)
 new <- orsf_fit(x[train, ],
                 y[train, ],
                 n_tree = 500,
                 mtry_ = 3,
                 leaf_min_events_ = 3,
                 leaf_min_obs_ = 15,
                 oobag_pred_ = TRUE,
                 cph_iter_max_ = 10,
                 cph_pval_max_ = .8)

 print(Score(list(new = 1-new$surv_oobag),
             formula = Surv(time, status) ~ 1,
             data = data.frame(as.matrix(y[train,])),
             times = 1000,
             summary = 'IPA'))

 rfs <- rfsrc(Surv(time, status) ~ .,
              nodesize = 15,
              data = as.data.frame(cbind(y[train,], x[train,])))

 pred_new <- orsf_pred_uni(new$forest, x[-train,], time_dbl = 1000)
 pred_rfs <- predictRisk(rfs, newdata = as.data.frame(x[-train,]), times = 1000)

 sc <- Score(
  object = list(new = pred_new,
                old = pred_rfs),
  formula = Surv(time, status) ~ 1,
  data = data.frame(as.matrix(y)[-train,]),
  times = 1000,
  summary = 'IPA'
 )

 res[[i]] <- sc

}

library(tidyverse)

map(res, 'AUC') |>
 map('score') |>
 bind_rows() |>
 group_by(model) |>
 summarize(auc=mean(AUC))

map(res, 'Brier') |>
 map('score') |>
 bind_rows() |>
 group_by(model) |>
 summarize(bri=mean(Brier),
           ipa=mean(IPA))

res[[1]]$AUC$score

set.seed(2)
old <- orsf_fit(x, y, n_tree = 100, mtry_ = 2, leaf_min_events_ = 20, oobag_pred_ = TRUE)

library(riskRegression)

sc <- Score(
 object = list(new = 1-new$surv_oobag,
               old = 1-old$surv_oobag),formula = Surv(time, status) ~ 1,
 data = .pbc,
 times = 1000,
 summary = 'IPA'
)





all(new$forest[[1]]$leaf_nodes == old$forest[[1]]$leaf_nodes)
all(new$forest[[1]]$betas == old$forest[[1]]$betas)

pred_new <- orsf_pred_uni(new$forest, x, time_dbl = 1000)
pred_old <- orsf_pred_uni(old$forest, x, time_dbl = 1000)

all(pred_new == pred_old)
