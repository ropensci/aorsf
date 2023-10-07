
set.seed(329)

#' @srrstats {G5.0} *use pbc/flchain data, standards for survival analysis*

library(survival)

# flchain ----

data("flchain", package = 'survival')

flc <- flchain
flc$chapter <- NULL
flc <- na.omit(flc)
flc <- flc[flc$futime > 0, ]

names(flc)[names(flc) == 'futime'] <- 'time'
names(flc)[names(flc) == 'death'] <- 'status'

# make sorted x and y matrices for testing internal cpp functions
flc_mats <- prep_test_matrices(flc, outcomes = c("time", "status"))

# lung ----

lcd <- survival::lung
lcd$inst <- NULL
lcd <- na.omit(lcd)

# make sorted x and y matrices for testing internal cpp functions
lcd_mats <- prep_test_matrices(lcd, outcomes = c("time", "status"))

# pbc ----

pbc <- pbc_orsf

pbc$id <- NULL

pbc_status_12 <- pbc
pbc_status_12$status <- pbc_status_12$status + 1

pbc_scale <- pbc_noise <- pbc

vars <- c('bili', 'chol', 'albumin', 'copper', 'alk.phos', 'ast')

for(i in vars){
 pbc_noise[[i]] <- add_noise(pbc_noise[[i]])
 pbc_scale[[i]] <- change_scale(pbc_scale[[i]])
}

# make sorted x and y matrices for testing internal cpp functions
pbc_mats <- prep_test_matrices(pbc, outcomes = c("time", "status"))

# data lists ----

data_list_pbc <- list(pbc_standard = pbc,
                      pbc_status_12 = pbc_status_12,
                      pbc_scaled = pbc_scale,
                      pbc_noised = pbc_noise)

# matric lists ----

mat_list_surv <- list(pbc = pbc_mats,
                      flc = flc_mats,
                      lcd = lcd_mats)

# standard fits ----

# standards used to check validity of other fits

seeds_standard <- 329
n_tree_test <- 10

controls <- list(
 fast = orsf_control_fast(),
 cph = orsf_control_cph(),
 net = orsf_control_net(),
 custom = orsf_control_custom(beta_fun = f_pca)
)

fit_standard_pbc <- lapply(
 controls,
 function(cntrl){
  orsf(pbc,
       formula = time + status ~ .,
       n_tree = n_tree_test,
       control = cntrl,
       tree_seed = seeds_standard)
 }
)

# training and testing data ----

pred_types_surv <- c(risk = 'risk',
                     surv = 'surv',
                     chf = 'chf',
                     mort = 'mort',
                     leaf = 'leaf')

pbc_train_rows <- sample(nrow(pbc_orsf), size = 170)

pbc_train <- pbc[pbc_train_rows, ]
pbc_test <- pbc[-pbc_train_rows, ]


