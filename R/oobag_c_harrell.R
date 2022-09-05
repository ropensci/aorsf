

#' Harrell's C-statistic
#'
#' This function is for testing and internal use.
#'
#' @param y_mat outcome matrix
#' @param s_vec vector of predicted survival
#'
#' @return the C-statistic
#'
#' @noRd
#'

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
