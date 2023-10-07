


# flc %>% dim()
#
# microbenchmark::microbenchmark(
#  aorsf = orsf(flc, time + status ~ .,
#               n_thread = 0,
#               n_split = 10,
#               leaf_min_obs = 10),
#  rfsrc = randomForestSRC::rfsrc(Surv(time, status) ~ .,
#                                 data = flc,
#                                 nodesize = 10,
#                                 samptype = 'swr'),
#  times = 1
# )
#
# microbenchmark::microbenchmark(
#  aorsf = orsf(pbc,
#               time + status ~ .,
#               n_thread = 0,
#               n_split = 10,
#               leaf_min_obs = 10),
#  rfsrc = randomForestSRC::rfsrc(Surv(time, status) ~ .,
#                                 data = pbc,
#                                 nodesize = 10,
#                                 samptype = 'swr'),
#  ranger = ranger::ranger(Surv(time, status) ~ .,
#                          data = pbc,
#                          num.thread = 0,
#                          min.node.size = 10),
#  times = 10
# )
