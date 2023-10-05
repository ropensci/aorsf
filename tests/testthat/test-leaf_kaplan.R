
#' @srrstats {G5.4} **Correctness tests** *test that statistical algorithms produce expected results to some fixed test data sets. I use the flchain data and compare the aorsf kaplan meier routine to that of the survival package.*

#' @srrstats {G5.4b} *Correctness tests include tests against previous implementations, explicitly calling those implementations in testing.*

#' @srrstats {G5.5} *Correctness tests are run with a fixed random seed*
set.seed(329)

weights <- sample(1:5, nrow(flc), replace = TRUE)

# fit a normal tree with no bootstrap weights
fit <- orsf(flc,
            futime + death ~ .,
            n_tree = 1,
            weights = weights,
            tree_seeds = 1,
            oobag_pred_type = 'none',
            # this makes every observation part of the in-bag data
            sample_fraction = 1,
            sample_with_replacement = FALSE,
            # this forces the tree to make a leaf at the root
            split_rule = 'cstat',
            split_min_stat = 0.999)

# so the result should be equivalent to fitting a kaplan-meier curve
# to the original training data, using replicate weights
aorsf_surv <- fit$forest$leaf_pred_prob[[1]][[1]]
aorsf_time <- fit$forest$leaf_pred_indx[[1]][[1]]

kap <- survfit(Surv(futime, death) ~ 1, data = flc, weights = weights)

test_that(
 desc = 'aorsf kaplan has same time values as survfit',
 code = {expect_equal(kap$time, aorsf_time, tolerance = 1e-9)}
)

test_that(
 desc = 'aorsf kaplan has same surv values as survfit',
 code = {expect_equal(kap$surv, aorsf_surv, tolerance = 1e-9)}
)
