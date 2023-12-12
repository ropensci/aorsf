# verbosity prints grow, predict, and importance notes

    Code
      fit_verbose <- orsf(pbc, time + status ~ ., verbose_progress = TRUE, n_tree = n_tree_test,
      importance = "negate")
    Output
      Growing trees: 100%. 
      Computing predictions: 100%. 
      Computing importance: 100%. 

---

    Code
      fit_verbose <- orsf(pbc, time + status ~ ., verbose_progress = TRUE, n_tree = n_tree_test,
      importance = "negate", n_thread = 1)
    Output
      Growing trees: 100%. 
      Computing predictions: 100%. 
      Computing importance: 100%. 

