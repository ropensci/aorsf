
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

# penguins ----

penguins <- penguins_orsf

penguins_binary <- penguins

penguins_binary$species <- factor(
 penguins_binary$species,
 levels = c("Adelie", "Chinstrap", "Gentoo"),
 labels = c("Adelie", "not_Adelie", "not_Adelie")
)

penguins_scale <- penguins_noise <- penguins

vars <- c("bill_length_mm", "bill_depth_mm",
          "flipper_length_mm", "body_mass_g")

for(i in vars){
 penguins_noise[[i]] <- add_noise(penguins_noise[[i]])
 penguins_scale[[i]] <- change_scale(penguins_scale[[i]])
}

# make sorted x and y matrices for testing internal cpp functions
penguins_mats <- prep_test_matrices(penguins, outcomes = c("species"))

# data lists ----

data_list_pbc <- list(pbc_standard = pbc,
                      pbc_status_12 = pbc_status_12,
                      pbc_scaled = pbc_scale,
                      pbc_noised = pbc_noise)

data_list_penguins <- list(penguins_standard = penguins,
                           penguins_binary = penguins_binary,
                           penguins_scaled = penguins_scale,
                           penguins_noised = penguins_noise)


# matrix lists ----

mat_list_surv <- list(pbc = pbc_mats,
                      flc = flc_mats,
                      lcd = lcd_mats)

# standard fits ----

# standards used to check validity of other fits

seeds_standard <- 329
n_tree_test <- 5

controls_surv <- list(
 fast = orsf_control_survival(method = 'glm', scale_x = FALSE, max_iter = 1),
 net = orsf_control_survival(method = 'net'),
 custom = orsf_control_survival(method = f_pca)
)

fit_standard_pbc <- lapply(
 controls_surv,
 function(cntrl){
  orsf(pbc,
       formula = time + status ~ .,
       n_tree = n_tree_test,
       control = cntrl,
       tree_seed = seeds_standard)
 }
)

controls_clsf <- list(
 fast = orsf_control_classification(method = 'glm', scale_x = FALSE, max_iter = 1),
 net = orsf_control_classification(method = 'net'),
 custom = orsf_control_classification(method = f_pca)
)

fit_standard_penguin_species <- lapply(
 controls_clsf,
 function(cntrl){
  orsf(penguins,
       formula = species ~ .,
       n_tree = n_tree_test,
       control = cntrl,
       tree_seed = seeds_standard)
 }
)

controls_regr <- list(
 fast = orsf_control_regression(method = 'glm', scale_x = FALSE, max_iter = 1),
 net = orsf_control_regression(method = 'net'),
 custom = orsf_control_regression(method = f_pca)
)

fit_standard_penguin_bills <- lapply(
 controls_regr,
 function(cntrl){
  orsf(penguins,
       formula = bill_length_mm ~ .,
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
                     leaf = 'leaf',
                     time = 'time')

pred_types_clsf <- c(prob = 'prob',
                     class = 'class',
                     leaf = 'leaf')

pred_types_regr <- c(mean = 'mean',
                     leaf = 'leaf')

pbc_train_rows <- sample(nrow(pbc_orsf), size = 170)
pbc_train <- pbc[pbc_train_rows, ]
pbc_test <- pbc[-pbc_train_rows, ]

penguins_train_rows <- sample(nrow(penguins_orsf), size = 180)
penguins_train <- penguins[penguins_train_rows, ]
penguins_test <- penguins[-penguins_train_rows, ]

penguins_binary_train <- penguins_binary[penguins_train_rows, ]
penguins_binary_test <- penguins_binary[-penguins_train_rows, ]

