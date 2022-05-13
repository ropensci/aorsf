

oobag_c_harrell <- function(y_mat, s_vec){

 sorted <- order(y_mat[, 1], -y_mat[, 2])

 y_mat <- y_mat[sorted, ]
 s_vec <- s_vec[sorted]

 time = y_mat[, 1]
 status = y_mat[, 2]
 events = which(status == 1)

 k = nrow(y_mat)

 total <- 0
 concordant <- 0

 for(i in events){

  if(i+1 <= k){

   for(j in seq(i+1, k)){

    if(time[j] > time[i]){

     total <- total + 1

     if(s_vec[j] > s_vec[i]){

      concordant <- concordant + 1

     } else if (s_vec[j] == s_vec[i]){

      concordant <- concordant + 0.5

     }

    }

   }

  }

 }

 concordant / total

}

fit_custom_oobag <- orsf(pbc_orsf,
                         formula = Surv(time, status) ~ . - id,
                         oobag_fun = oobag_c_harrell,
                         n_tree = 10,
                         tree_seeds = 1:10)

fit_standard_oobag <- orsf(pbc_orsf,
                           formula = Surv(time, status) ~ . - id,
                           n_tree = 10,
                           tree_seeds = 1:10)

test_that(
 desc = 'tree seeds show that a custom oobag fun matches the internal one',
 code = {
  expect_equal(
   fit_standard_oobag$eval_oobag$stat_values,
   fit_custom_oobag$eval_oobag$stat_values
  )
 }
)



