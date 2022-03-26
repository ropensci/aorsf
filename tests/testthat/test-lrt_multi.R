
#' @srrstatsTODO {G5.4} **Correctness tests** *test that statistical algorithms produce expected results to some fixed test data sets. I simulate arbitrary data and compare the aorsf likelihood ratio test to the same algorithm used in survival::survdiff().*
#'
set.seed(329)

leaf_min_events <- 1
leaf_min_obs    <- 15
n_split         <- 5
n_total         <- 100

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

lrt_multi_vals <- lapply(
 X = list(ctns = XB_ctns,
          catg = XB_catg,
          bnry = XB_bnry),
 FUN = function(XB){
  lrt_multi_testthat(y_node_ = y,
                     w_node_ = w,
                     XB_ = XB,
                     n_split_ = n_split,
                     leaf_min_events_ = leaf_min_events,
                     leaf_min_obs_ = leaf_min_obs)
 }
)

test_values = lrt_multi_vals[[1]]

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

 cp_first = xb_uni[min(which(cp_stats$valid_cp))]
 cp_last  = xb_uni[max(which(cp_stats$valid_cp))]

 index_last <- max(which(test_values$cutpoints!=0))

 test_that(
  desc = 'same chi-squared stats as survival survdiff',
  code = {

   expect_equal(cp_first, test_values$cutpoints[1])
   expect_equal(cp_last, test_values$cutpoints[index_last])

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

run_lrt_multi_tests(lrt_multi_vals$ctns, XB_ctns)
run_lrt_multi_tests(lrt_multi_vals$catg, XB_catg)
run_lrt_multi_tests(lrt_multi_vals$bnry, XB_bnry)


