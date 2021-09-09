

# x is binary, with 3 cases
# -- if leaf_min_obs is 2, left node has enough but right doesn't
x_bnry <- c(0, 1, 0, 1, 0, 0, 0, 0, 0, 1)

# almost all observations had an event
y <- matrix(
 c(1,2,3,4,5,6,7,8,9,10,
   0,1,1,1,1,1,1,1,1,1),
 byrow = FALSE,
 ncol = 2
)

colnames(y) <- c("time", "status")

cp <- rep(Inf, 3)

# if 0 is the cut-point,
# -- 6 events and 6 obs in the left node
# -- 3 events and 3 obs in the right now
# -- valid cut-point with leaf_min_obs <= 2 and leaf_min_events <= 2
# if 1 is the cut-point, everything is in the left node
# -- invalid cut-point

find_cutpoints(
   cp,
   lc = x_bnry,
   y = y,
   leaf_min_obs = 1,
   leaf_min_events = 1
)

test_that(desc = 'binary first case is correct',
          code = {expect_equal(cp, c(0, Inf, Inf))})

find_cutpoints(
 cp,
 lc = x_bnry,
 y = y,
 leaf_min_obs = 2,
 leaf_min_events = 2
)


test_that(desc = 'binary second case is correct',
          code = {expect_equal(cp, c(0, Inf, Inf))})

# still works when we have exactly 3 and 3 events and obs required
find_cutpoints(
   cp,
   lc = x_bnry,
   y = y,
   leaf_min_obs = 3,
   leaf_min_events = 3
)

test_that(desc = 'binary third case is correct',
          code = {expect_equal(cp, c(0, Inf, Inf))})

# no valid cutpoints if leaf_min_obs is 4
find_cutpoints(
   cp,
   lc = x_bnry,
   y = y,
   leaf_min_obs = 4,
   leaf_min_events = 3
)

test_that(desc = 'binary fourth case is correct',
          code = {expect_equal(cp, c(Inf, Inf, Inf))})

# no valid cutpoints if leaf_min_obs is 4
find_cutpoints(
   cp,
   lc = x_bnry,
   y = y,
   leaf_min_obs = 3,
   leaf_min_events = 4
)

test_that(desc = 'binary fifth case is correct',
          code = {expect_equal(cp, c(Inf, Inf, Inf))})

# x is categorical, with 3 categories
# -- if leaf_min_obs is 2, left node has enough but right doesn't
x_catg <- c(0, 1, 2, 1, 2, 0, 2, 0, 2, 1)


# 2 events, 3 obs in 0 catg
# 3 events, 3 obs in 1 catg
# 4 events, 4 obs in 2 catg
# --> 0 and 1 are valid cp's if leaf_min_obs = 2 and leaf_min_event = 2
# --> 0 and 1 are valid cp's if leaf_min_obs = 3 and leaf_min_event = 2
# --> 1 is the only valid cp if leaf_min_obs = 3 and leaf_min_event = 3
#
find_cutpoints(
   cp,
   lc = x_catg,
   y = y,
   leaf_min_obs = 1,
   leaf_min_events = 1
)

test_that(desc = 'categorical first case is correct',
          code = {expect_equal(cp, c(0, 1, Inf))})

find_cutpoints(
   cp,
   lc = x_catg,
   y = y,
   leaf_min_obs = 2,
   leaf_min_events = 2
)

test_that(desc = 'categorical second case is correct',
          code = {expect_equal(cp, c(0, 1, Inf))})

find_cutpoints(
   cp,
   lc = x_catg,
   y = y,
   leaf_min_obs = 3,
   leaf_min_events = 2
)

test_that(desc = 'categorical third case is correct',
          code = {expect_equal(cp, c(0, 1, Inf))})

find_cutpoints(
   cp,
   lc = x_catg,
   y = y,
   leaf_min_obs = 3,
   leaf_min_events = 3
)

test_that(desc = 'categorical fourth case is correct',
          code = {expect_equal(cp, c(Inf, 1, Inf))})




# x is continuous
x_ctns <- c(10, 9, 9, 8, 7, 5, 6, 3, 2, 1)

y <- matrix(
   c(1,2,3,4,5,6,7,8,9,10,
     0,1,1,1,1,1,1,1,1,1),
   byrow = FALSE,
   ncol = 2
)

cp <- rep(0, 2)

find_cutpoints(
   cp,
   lc = x_ctns,
   y = y,
   leaf_min_obs = 2,
   leaf_min_events = 2
)

weights <- rep(1, nrow(y))

find_cutpoints_ctns(cp, x_ctns, y, weights, 2, 2)

# # benchmark
#
# fit <- survival::coxph(flchain_y ~ flchain_x[, c('age', 'kappa')])
# x <- predict(fit)
# cp <- rep(Inf, 5)
#
# # continuous
# microbenchmark::microbenchmark(
#  c = find_cutpoints(cp, x, flchain_y, leaf_min_obs = 3, leaf_min_events = 6),
#  r = quantile(x, probs = c(0.05, 0.25, 0.50, 0.75, 0.95))
# )
#
# fit <- survival::coxph(flchain_y ~ flchain_x[, c('sexF')])
# x <- predict(fit)
# cp <- rep(Inf, 5)
#
# # binary
# microbenchmark::microbenchmark(
#    c = find_cutpoints(cp, x, flchain_y, leaf_min_obs = 3, leaf_min_events = 6),
#    r = quantile(x, probs = c(0.05, 0.25, 0.50, 0.75, 0.95))
# )
#
# fit <- survival::coxph(flchain_y ~ flchain_x[, c('sexF', 'mgus')])
# x <- predict(fit)
# cp <- rep(Inf, 5)
#
# # categorical (4 categories)
# microbenchmark::microbenchmark(
#    c = find_cutpoints(cp, x, flchain_y, leaf_min_obs = 3, leaf_min_events = 6),
#    r = quantile(x, probs = c(0.05, 0.25, 0.50, 0.75, 0.95))
# )





