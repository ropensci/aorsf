

orsf_fit_xy <- function(x,
                        y,
                        ntree = 1000,
                        nsplit = 5,
                        mtry = 4,
                        leaf_min_events = 5,
                        leaf_min_obs = 10,
                        cph_method = 1,
                        cph_eps = 1e-2,
                        cph_iter_max = 5,
                        cph_prune_thresh = 0.15,
                        cph_rescale = TRUE){

 # data must be sorted from low to high by time
 sorted_by_time <- frank(y[, 1], ties.method = 'random')

 orsf_fit(x[sorted_by_time, ],
          y[sorted_by_time, ],
          ntree,
          mtry,
          nsplit,
          leaf_min_events,
          leaf_min_obs,
          cph_method,
          cph_eps,
          cph_iter_max,
          cph_prune_thresh,
          cph_rescale)

}
