
library(survival)
devtools::load_all()

.pbc <-  pbc[order(pbc$time), ]
.pbc <- .pbc[complete.cases(.pbc), ]
.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1

x <- as.matrix(.pbc[, c('trt','age','ascites','edema','bili')])
y <- Surv(.pbc$time, .pbc$status)

# x <- flchain_x
# y <- flchain_y

set.seed(1)

tree = ostree_fit_arma(
  x,
  y,
  mtry = 4,
  n_cps = 5,
  leaf_min_events = 10,
  leaf_min_obs = 10,
  verbose = T
)

ne
