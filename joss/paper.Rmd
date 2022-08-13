---
# Example from https://joss.readthedocs.io/en/latest/submitting.html
title: 'aorsf: An R package for supervised learning using the oblique random survival forest'
tags:
  - R
  - machine learning
  - supervised learning
  - survival
  - random forest
authors:
  - name: Byron C. Jaeger
    orcid: 0000-0001-7399-2299
    affiliation: 1 # (Note: multiple affiliations must be quoted)
  - name: Sawyer Welden
    affiliation: 1
  - name: Kristin Lenoir
    affiliation: 1
  - name: Nicholas M Pajewski
    affiliation: 1
affiliations:
 - name: Wake Forest University School of Medicine
   index: 1
citation_author: Jaeger et. al.
date: 11 May 2022
year: 2022
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

# Summary

The random survival forest (RSF) is a supervised learning method for right-censored time-to-event data that combines predictions from a large set of survival decision trees [@ishwaran2008random]. Similar to random forests for classification and regression [@breiman2001random], trees in the RSF are de-correlated by using random subsets of training data to grow each tree and considering a random subset of predictor variables at each non-terminal node. Decision trees in the RSF may be axis based or oblique. Axis-based trees split non-terminal nodes using individual predictor variables, whereas oblique trees use a linear combination of variables. Although using oblique instead of axis based decision trees to grow the RSF improves its Brier score and concordance index in prediction tasks, the computational overhead of fitting an oblique RSF may be substantially higher than an axis based RSF [@jaeger2019oblique], making it difficult to use the oblique RSF in applied settings. 


``aorsf`` is an R package designed to develop and compute predictions with oblique RSFs efficiently. Extensions of ``aorsf``'s core features are supported by allowing users to supply their own function to identify a linear combination of inputs when growing oblique survival trees. The target audience includes both __analysts__ aiming to develop an accurate risk prediction model (e.g., see @segar2021development) and __researchers__ who want to conduct experiments comparing different techniques for identifying linear combinations of predictor variables (e.g., see @katuwal2020heterogeneous). Key features of ``aorsf`` include computational efficiency compared to existing software, extensive unit and integration testing that ensure cross-platform consistency and reproducibility, and user-friendly documentation paired with an application programming interface that facilitates proper usage of the core algorithms.

# Existing software 

The `obliqueRF` and `RLT` R packages support oblique classification and regression random forests but not oblique RSFs. The `ranger` and `randomForestSRC` packages support axis based RSFs but not oblique RSFs. The ``obliqueRSF`` R package fits oblique RSFs, but its computational inefficiency is a barrier in applied settings. ``aorsf`` is designed for computational efficiency and flexibility, allowing users to grow oblique RSFs using Newton Raphson scoring (the default), penalized Cox regression (the method used by `obliqueRSF`), or a user-defined function to identify linear combinations of variables. 

# Newton Raphson scoring

The default routine for creating linear combinations of predictor variables in ``aorsf`` applies Newton Raphson scoring to the partial likelihood function of the Cox regression model. ``aorsf`` uses the same approach as the `survival` package to complete this estimation procedure efficiently. Full details on the steps involved have been made available by @therneau_survival_2022. Briefly, a vector of estimated regression coefficients, $\hat{\beta}$, is updated in each step of the procedure based on its first derivative, $U(\hat{\beta})$, and second derivative, $H(\hat{\beta})$: 

$$ \hat{\beta}^{k+1} =  \hat{\beta}^{k} + U(\hat{\beta} = \hat{\beta}^{k})\, H^{-1}(\hat{\beta} = \hat{\beta}^{k})$$
While it is standard practice in statistical modeling to iteratively update $\hat{\beta}$ until a convergence threshold is met, the default approach in ``aorsf`` only completes one iteration (see `?orsf_control_cph`). The rationale for this approach is based on two points. First, while completing more iterations reduces bias in the regression coefficients, general benchmarking experiments have found it does not improve the discrimination or Brier score of the oblique RSF [@aorsf_arxiv]. Second, computing $U$ and $H$ requires computation and exponentiation of the vector $X\hat{\beta}$, where $X$ is the matrix of predictor values, but these steps can be skipped on the first iteration if an initial value of $\hat{\beta} = 0$ is chosen, allowing for a reduction in required computation.

# Benchmarking 

The increased efficiency of ``aorsf`` versus `obliqueRSF` results from improved memory management and Newton Raphson scoring. In general benchmarks of prediction accuracy and computational efficiency, ``aorsf`` has matched or exceeded the prediction accuracy of `obliqueRSF` while running two orders of magnitude faster [@aorsf_arxiv]. 

# Acknowledgements

The development of this software was supported by the Center for Biomedical Informatics at Wake Forest University School of Medicine and the National Center for Advancing Translational Sciences (NCATS), National Institutes of Health (NIH), through Grant Award Number UL1TR001420. Dr. Pajewski was additionally supported by grant award P30 AG021332 from the NIH. The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH

# References
