
#' @srrstats {G5.4} **Correctness tests** *test that statistical algorithms produce expected results to some fixed test data sets. I simulate arbitrary data and compare the aorsf likelihood ratio test to the same algorithm used in survival::survdiff().*

#' @srrstats {G5.4b} *Correctness tests include tests against previous implementations, explicitly calling those implementations in testing.*

#' @srrstats {G5.5} *Correctness tests are run with a fixed random seed*
set.seed(329555)

#' @srrstats {G5.6} **Parameter recovery tests** *the likelihood ratio test returns expected values consistent with the survival implementation for randomly generated data*

#' @srrstats {ML7.8} *Explicitly test my implementation of the likelihood ratio test, used as a loss function when determining the best split for a given node. I do not test other loss functions because this is the only loss function that aorsf implements.*

run_lrt_multi_tests <- function(test_values, XB){

 xb_uni <- unique(XB)

 cp_stats <-
  sapply(
   X = xb_uni,
   FUN = function(x){
    c(
     cp = x,
     e_right = sum(status[XB > x]),
     e_left = sum(status[XB <= x]),
     n_right = sum(XB > x),
     n_left = sum(XB <= x)
    )
   }
  )

 cp_stats <- as.data.frame(t(cp_stats))

 cp_stats$valid_cp = with(
  cp_stats,
  e_right >= leaf_min_events & e_left >= leaf_min_events  &
   n_right >= leaf_min_obs & n_left >= leaf_min_obs
 )

 if(!any(cp_stats$valid_cp)){
  return(NULL)
 }

 cp_first = xb_uni[min(which(cp_stats$valid_cp))]
 cp_last  = xb_uni[max(which(cp_stats$valid_cp))]

 test_that(
  desc = 'same chi-squared stats as survival survdiff',
  code = {

   for(i in seq(index_last)){

    XB_cut <- as.numeric(XB > test_values$cutpoints[i])

    chisq <-
     getElement(survival::survdiff(survival::Surv(time, status) ~ XB_cut),
                'chisq')

    # same chi square stat
    expect_equal(test_values$statistic[i], chisq)

    # valid event and observation counts
    e_right <- sum(status[XB > test_values$cutpoints[i]])
    n_right <- sum(XB > test_values$cutpoints[i])

    e_left <- sum(status[XB <= test_values$cutpoints[i]])
    n_left <- sum(XB <= test_values$cutpoints[i])

    expect_equal(e_left + e_right, n_event)
    expect_equal(n_left + n_right, n_total)

    expect_true(e_right >= leaf_min_events)
    expect_true(e_left >= leaf_min_events)

    expect_true(n_right >= leaf_min_obs)
    expect_true(n_left >= leaf_min_obs)

   }

  }

 )

}


n_total         <- 250

.leaf_min_events <- c(1, 3, 5, 10, 15)

for(leaf_min_events in .leaf_min_events){

 leaf_min_obs    <- leaf_min_events + 10

 XB_ctns <- sort(rnorm(n_total))
 XB_catg <- round(XB_ctns)
 XB_bnry <- as.numeric(XB_ctns > 0)

 XB_catg <- XB_catg + abs(min(XB_catg)) + 1
 XB_bnry <- XB_bnry + 1

 prob <- (XB_ctns + abs(min(XB_ctns))) / (max(XB_ctns)+abs(min(XB_ctns)))

 status  <- rbinom(n = n_total, prob = prob, size = 1)

 n_event <- sum(status)

 time    <- seq(n_total, 1)
 t_sort  <- order(time)
 status  <- status[t_sort]
 XB_ctns <- XB_ctns[t_sort]
 XB_catg <- XB_catg[t_sort]
 XB_bnry <- XB_bnry[t_sort]
 time    <- time[t_sort]

 y <- cbind(time=time, status=status)
 w <- rep(1, n_total)

 cp_bounds <- lapply(
  X = list(ctns = XB_ctns,
           catg = XB_catg,
           bnry = XB_bnry),
  FUN = function(XB){
   cp_find_bounds_R(y_node = y,
                    w_node = w,
                    XB = XB,
                    xb_uni = unique(XB),
                    leaf_min_events = leaf_min_events,
                    leaf_min_obs = leaf_min_obs)
  }
 ) %>%
  lapply(subset, valid_cp)

 for(i in seq_along(cp_bounds)){

  XB <- switch (names(cp_bounds)[i],
                'ctns' = XB_ctns,
                'catg' = XB_catg,
                'bnry' = XB_bnry)

  for(j in seq_along(cp_bounds[[i]]$cp)){


   group <- as.numeric(XB > cp_bounds[[i]]$cp[j])


   expect_equal(
    node_compute_lrt_exported(y, w, group),
    survival::survdiff(survival::Surv(time, status) ~ group)$chisq
   )

  }


 }

}


# # benchmark does not need to be tested every time
#
# bm <- microbenchmark::microbenchmark(
#  R = survival::survdiff(survival::Surv(time, status) ~ group)$chisq,
#  cpp = node_compute_lrt_exported(y, w, group),
#  times = 50
# )
#
# expect_lt(
#  median(bm$time[bm$expr == 'cpp']),
#  median(bm$time[bm$expr == 'R'])
# )



