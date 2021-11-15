
library(survival)

.pbc <- pbc[complete.cases(pbc), ]
.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1
.pbc$stage <- factor(.pbc$stage, ordered = TRUE)

pbc_orsf <- .pbc

usethis::use_data(pbc_orsf, overwrite = TRUE)
