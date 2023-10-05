
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

# data lists ----

data_list_pbc <- list(pbc_standard = pbc,
                      pbc_status_12 = pbc_status_12,
                      pbc_scaled = pbc_scale,
                      pbc_noiced = pbc_noise)

# standards used to check validity of other fits

# standard fits ----

seeds_standard <- c(5, 20, 1000, 30, 50, 98, 22, 100, 329, 10)

fit_standard_pbc <- orsf(pbc,
                         formula = time + status ~ .,
                         n_tree = 10,
                         tree_seed = seeds_standard)
