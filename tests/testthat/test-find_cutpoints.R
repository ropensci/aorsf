
test_that(
 desc = 'cutpoints are unique and correct',
 code = {

  skip_on_cran()

  for(i in seq_along(mat_list_surv)){

   y <- mat_list_surv[[i]]$y
   w <- mat_list_surv[[i]]$w

   if(nrow(y) > 100) {y <- y[1:100, ]; w <- w[1:100]}

   for(cp_type in c("ctns", "bnry", "catg")){

    xb <- switch(
     cp_type,
     'ctns' = rnorm(nrow(y)),
     'bnry' = rbinom(nrow(y), size = 1, prob = 1/2),
     'catg' = rbinom(nrow(y), size = 5, prob = 1/2)
    )

    xb_uni <- unique(xb)
    # leaf_min_events <- 5
    # leaf_min_obs <- 10
    for(leaf_min_events in c(1, 5)){

     for(leaf_min_obs in c(leaf_min_events + c(0, 5))){

      cp_stats <- cp_find_bounds_R(y, w, xb, xb_uni, leaf_min_events, leaf_min_obs)

      cp_index <- find_cuts_survival_exported(y, w, xb,
                                              leaf_min_events,
                                              leaf_min_obs,
                                              split_rule_R = 1)


      cps_r <- cp_stats$cp[cp_stats$valid_cp]

      cps_cpp <- sort(xb)[cp_index$cuts_all+1]

      expect_equal(length(cps_cpp), length(unique(cps_cpp)))
      expect_true(is_equivalent(cps_r, cps_cpp))
      expect_true(all(cp_index$cuts_sampled %in% cp_index$cuts_all))
      expect_equal(unique(cp_index$cuts_sampled), cp_index$cuts_sampled)

      cps_sampled <- sort(xb)[cp_index$cuts_sampled+1]

      g_list <- lapply(
       cps_sampled, function(cp){as.numeric(xb <= cp)}
      )

      logrank_stats <- sapply(
       g_list, function(gg){compute_logrank_exported(y, w, gg)}
      )

      expect_equal(cps_sampled[which.max(logrank_stats)], cp_index$best_cut)

     }

    }

   }

  }

 }
)




# benchmark does not need to be tested every time

# bm <- microbenchmark::microbenchmark(
#
#   cp_stats = cp_find_bounds_R(y, w, xb, xb_uni, leaf_min_events, leaf_min_obs),
#
#   cp_index = find_cutpoints_survival_exported(y, w, xb,
#                                                leaf_min_events,
#                                                leaf_min_obs),
#
#  times = 50
#
# )
#
# expect_lt(
#  median(bm$time[bm$expr == 'cpp']),
#  median(bm$time[bm$expr == 'R'])
# )




