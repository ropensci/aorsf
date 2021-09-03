
library(survival)

.pbc <-  pbc[order(pbc$time), ]
.pbc <- .pbc[complete.cases(.pbc), ]
.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1

x <- as.matrix(.pbc[, c('trt','age','ascites','edema','bili')])
y <- Surv(.pbc$time, .pbc$status)

set.seed(1)

tree = ostree_fit_arma(
 x,
 y,
 mtry = 3,
 n_vars_lc = 2,
 n_cps = 3,
 leaf_min_events = 5,
 leaf_min_obs = 10,
 verbose = TRUE
)

tmp = ostree_predict(tree, x)

tmp < tree$node_0$cutpoint
