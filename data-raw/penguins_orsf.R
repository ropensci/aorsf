library(palmerpenguins)

penguins_orsf <- stats::na.omit(penguins) |>
 dplyr::mutate(bill_length_mm = as.numeric(bill_length_mm),
               flipper_length_mm = as.numeric(flipper_length_mm))

usethis::use_data(penguins_orsf, overwrite = TRUE)
