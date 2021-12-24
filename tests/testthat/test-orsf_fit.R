

object = orsf(formula = Surv(time, status) ~ . - id,
              data_train = pbc_orsf,
              mtry = 5,
              n_split = 5,
              cph_do_scale = TRUE,
              n_tree = 500,
              leaf_min_obs = 10,
              importance = FALSE,
              oobag_pred = TRUE)

