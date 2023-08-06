

run_cp_bounds_test <- function(test_values, XB){

 xb_uni <- unique(XB)

 cp_stats <- cp_find_bounds_R(y, w, XB, xb_uni, leaf_min_events, leaf_min_obs)

 if(!any(cp_stats$valid_cp)){
  return(NULL)
 }

 cps_true_values <- sort(xb_uni[cp_stats$valid_cp])
 cps_test_values <- XB[test_values$XB_sorted+1][test_values$cp_index+1]

 test_that(
  desc = 'cutpoints identified are unique and valid',
  code = {

   expect_equal(
    length(cps_test_values), length(unique(cps_test_values))
   )

   expect_equal(cps_true_values, cps_test_values)

  }

 )

 test_that(
  desc = "group values are filled corresponding to the given cut-point",
  code = {

   group_cpp <- rep(0, length(XB))
   XB_sorted <- order(XB)-1

   for(i in seq_along(cps_true_values)){

    group_R = XB <= cps_true_values[i]

    if(i == 1) start <- 0 else start <- test_values$cp_index[i-1]+1

    node_fill_group_exported(
     group = group_cpp,
     XB_sorted = XB_sorted,
     start = start,
     stop = test_values$cp_index[i],
     value = 1
    )

    expect_equal(as.numeric(group_R),
                 as.numeric(group_cpp))

   }
  }
 )

}

.leaf_min_events <- c(1, 5, 50, nrow(pbc_orsf))

# leaf_min_events = 1

for(leaf_min_events in .leaf_min_events){

 leaf_min_obs    <- leaf_min_events + 10

 XB_ctns <- pbc_orsf$age
 XB_catg <- round(pbc_orsf$bili)
 XB_bnry <- as.numeric(pbc_orsf$sex)

 status  <- pbc_orsf$status
 time    <- pbc_orsf$time

 t_sort  <- order(time)
 status  <- status[t_sort]
 XB_ctns <- XB_ctns[t_sort]
 XB_catg <- XB_catg[t_sort]
 XB_bnry <- XB_bnry[t_sort]
 time    <- time[t_sort]

 y <- cbind(time=time, status=status)
 w <- rep(1, nrow(pbc_orsf))

 cp_bounds <- lapply(
  X = list(ctns = XB_ctns,
           catg = XB_catg,
           bnry = XB_bnry),
  FUN = function(XB){
   node_find_cps_exported(y_node = y,
                          w_node = w,
                          XB = XB,
                          leaf_min_events = leaf_min_events,
                          leaf_min_obs = leaf_min_obs)
  }
 )

 run_cp_bounds_test(test_values = cp_bounds$ctns, XB = XB_ctns)
 run_cp_bounds_test(cp_bounds$catg, XB = XB_catg)
 run_cp_bounds_test(cp_bounds$bnry, XB = XB_bnry)



}


# benchmark does not need to be tested every time

# bm <- microbenchmark::microbenchmark(
#
#  R = {
#   xb_uni = unique(XB_ctns)
#   cp_find_bounds_R(y_node = y,
#                    w_node = w,
#                    XB = XB_ctns,
#                    xb_uni = xb_uni,
#                    leaf_min_events = 5,
#                    leaf_min_obs = 10)
#  },
#
#  cpp = cp_find_bounds_exported(y_node = y,
#                                w_node = w,
#                                XB = XB_ctns,
#                                leaf_min_events = 5,
#                                leaf_min_obs = 10),
#
#  times = 50
#
# )
#
# expect_lt(
#  median(bm$time[bm$expr == 'cpp']),
#  median(bm$time[bm$expr == 'R'])
# )



