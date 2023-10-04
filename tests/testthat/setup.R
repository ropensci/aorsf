

# standard fit object used to check validity of other fits

seeds_standard <- c(5, 20, 1000, 30, 50, 98, 22, 100, 329, 10)

fit_standard <- orsf(pbc_orsf,
                     formula = time + status ~ . - id,
                     n_tree = 10,
                     tree_seed = seeds_standard)
