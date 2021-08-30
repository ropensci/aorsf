
library(survival)

df <- flchain

df$chapter <- NULL

time <- 'futime'
status <- 'death'

df_nomiss <- na.omit(df)

df_sorted <- df_nomiss[order(df_nomiss[[time]]),]

df_x <- df_sorted
df_x[[time]] <- NULL
df_x[[status]] <- NULL

flchain_x <- model.matrix(~.-1, data = df_x)

flchain_y <- Surv(time = df_sorted[[time]],
          event = df_sorted[[status]])

usethis::use_data(flchain_x, flchain_y, overwrite = TRUE)
