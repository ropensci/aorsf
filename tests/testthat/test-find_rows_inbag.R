
n_obs <- 100

rows_inbag_checker <- function(rows_oobag, n_obs){
 all_rows <- seq(n_obs)-1
 matrix(setdiff(all_rows, rows_oobag), ncol=1)
}

for(i in seq(100)){

 rows_oobag <- sort(sample(x = seq(0, n_obs-1), size = round(n_obs/3)))

 expect_equal(
  find_rows_inbag_exported(rows_oobag, n_obs),
  rows_inbag_checker(rows_oobag, n_obs)
 )

}
