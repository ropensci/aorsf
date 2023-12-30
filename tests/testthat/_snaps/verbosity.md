# verbosity prints grow, predict, and importance notes

    Code
      fit_verbose <- orsf(pbc, time + status ~ ., verbose_progress = TRUE, n_tree = n_tree_test,
      importance = "negate")

---

    Code
      fit_verbose <- orsf(pbc, time + status ~ ., verbose_progress = TRUE, n_tree = n_tree_test,
      importance = "negate", n_thread = 1)
    Output
      Growing trees: 100%. 
      Computing predictions: 100%. 
      Computing importance: 100%. 

# verbosity is carried by object

    Code
      pd <- orsf_pd_oob(fit_verbose, pred_spec_auto(island))

---

    Code
      vi <- orsf_vi(fit_verbose_2, importance = "negate")

---

    Code
      vs <- orsf_vs(fit_verbose)
    Output
      Selecting variables: 100%

