library(palmerpenguins)

penguins_orsf <- stats::na.omit(penguins)

usethis::use_data(penguins_orsf, overwrite = TRUE)
