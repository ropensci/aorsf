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

# verbosity is carried by object

    Code
      orsf_pd_oob(fit_verbose, pred_spec_auto(island))
    Output
      Computing dependence: 100%. 
             class    island      mean lwr        medn       upr
      1:    Adelie    Biscoe 0.4904928   0 0.470588235 1.0000000
      2:    Adelie     Dream 0.4169041   0 0.150000000 1.0000000
      3:    Adelie Torgersen 0.6145236   0 0.941176471 1.0000000
      4: Chinstrap    Biscoe 0.1191911   0 0.000000000 0.9666667
      5: Chinstrap     Dream 0.3739789   0 0.254562468 1.0000000
      6: Chinstrap Torgersen 0.2061961   0 0.007201859 1.0000000
      7:    Gentoo    Biscoe 0.3903161   0 0.051470588 1.0000000
      8:    Gentoo     Dream 0.2091170   0 0.000000000 0.9863014
      9:    Gentoo Torgersen 0.1792803   0 0.000000000 1.0000000

---

    Code
      orsf_vi(fit_verbose_2, importance = "negate")
    Output
      Computing importance: 100%. 
      flipper_length_mm       body_mass_g           species            island 
             0.52904478        0.34635451        0.31919457        0.14621939 
          bill_depth_mm               sex              year 
             0.13568507        0.09509735        0.01633449 

---

    Code
      orsf_vs(fit_verbose)
    Output
      Selecting variables: 100%
         n_predictors stat_value
      1:            3  0.9991496
      2:            4  0.9821797
      3:            5  0.9972835
      4:            6  0.9935440
      5:            7  0.9861669
                                                                                  predictors_included
      1:                                           island_Dream,bill_length_mm,bill_depth_mm,sex_male
      2:                               island_Dream,bill_length_mm,bill_depth_mm,body_mass_g,sex_male
      3:              island_Dream,island_Torgersen,bill_length_mm,bill_depth_mm,body_mass_g,sex_male
      4: island_Dream,island_Torgersen,bill_length_mm,bill_depth_mm,flipper_length_mm,body_mass_g,...
      5: island_Dream,island_Torgersen,bill_length_mm,bill_depth_mm,flipper_length_mm,body_mass_g,...
         predictor_dropped
      1:      island_Dream
      2:       body_mass_g
      3:  island_Torgersen
      4: flipper_length_mm
      5:              year

