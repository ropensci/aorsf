
library(survival)

.pbc <- pbc[complete.cases(pbc), ]
.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1
.pbc$stage <- factor(.pbc$stage, ordered = TRUE)

pbc_orsf <- .pbc

pbc_orsf$trt <- factor(pbc_orsf$trt,
                       levels = c(1, 2),
                       labels = c('d_penicill_main',
                                  'placebo'))
pbc_orsf$ascites <- factor(pbc_orsf$ascites)
pbc_orsf$hepato <- factor(pbc_orsf$hepato)
pbc_orsf$spiders <- factor(pbc_orsf$spiders)

usethis::use_data(pbc_orsf, overwrite = TRUE)
