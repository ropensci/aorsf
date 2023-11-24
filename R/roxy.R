# nocov start

# data allowed ------------------------------------------------------------

roxy_data_allowed <- function(){
 paste_collapse(
  x = c(
   "[data.frame]",
   "[tibble][tibble::tibble-package]",
   "[data.table][data.table::data.table-package]"
  ),
  sep = ', ',
  last = ' or '
 )
}

# Oblique Forest descriptor -----------------------------------------------

roxy_describe_ObliqueForest <- function(trained){
 paste("(*ObliqueForest*)",
       if(trained) "a trained" else "an",
       "oblique random forest object (see [orsf])")
}


# multi-threading ---------------------------------------------------------

roxy_n_thread_header <- function(action){
 paste("(_integer_) number of threads to use while ",
       action, ". Default is one thread. ",
       "To use the maximum number of threads that ",
       "your system provides for concurrent execution, ",
       "set `n_thread = 0`.", sep = "")
}

roxy_n_thread_details <- function(){
 "(_integer_) number of threads to use. Default is one thread."
}

# importance --------------------------------------------------------------

roxy_importance_header <- function(){
 "(_character_) Indicate method for variable importance:\n"
}

roxy_importance_none <- function(){
 "'none': no variable importance is computed."
}
roxy_importance_anova <- function(){
 "'anova': compute analysis of variance (ANOVA) importance"
}
roxy_importance_negate <- function(){
 "'negate': compute negation importance"
}
roxy_importance_permute <- function(){
 "'permute': compute permutation importance"
}

roxy_group_factors <- function(){
 "if `TRUE`, the importance of factor variables will be reported overall by aggregating the importance of individual levels of the factor. If `FALSE`, the importance of individual factor levels will be returned."
}

# oobag_fun ---------------------------------------------------------------

roxy_oobag_fun_header <- function(){

 "(_function_) to be used for evaluating out-of-bag prediction accuracy"

}

roxy_oobag_fun_default <- function(
  stat_label = c("- survival: Harrell's C-statistic (1982)",
                 "- classification: Area underneath the ROC curve (AUC-ROC)",
                 "- regression: Traditional prediction R-squared")
){
  paste("When `oobag_fun = NULL` (the default), the evaluation",
        "statistic is selected based on tree type\n\n",
        paste(stat_label, collapse = '\n'))
 }

roxy_oobag_fun_user <- function(){
 "if you use your own `oobag_fun` note the following:"
}

roxy_oobag_fun_inputs <- function(){
 "`oobag_fun` should have three inputs: `y_mat`, `w_vec`, and `s_vec`"
}

roxy_oobag_fun_ymat <- function(){
  paste("For survival trees, `y_mat` should be a two column matrix with",
  "first column named 'time' and second named 'status'. For classification",
  "trees, `y_mat` should be a matrix with number of columns = number of",
  "distinct classes in the outcome. For regression, `y_mat` should be a",
  "matrix with one column.")
}

roxy_oobag_fun_svec <- function(){
 "`s_vec` is a numeric vector containing predictions"
 # paste("`s_vec` is a numeric vector containing predictions. For survival,",
 #       "the number of columns should match the length of `oobag_pred_horison`.",
 #       "For classification, the number of columns should match the number of",
 #       "classes. For regression, the number of columns should be 1." )
}


roxy_oobag_fun_return <- function(){
 "`oobag_fun` should return a numeric output of length 1"
}

# oobag_fun ---------------------------------------------------------------

roxy_na_action_header <- function(data_label){

 paste0(
  "(_character_) what should happen when `", data_label,"` contains missing values (i.e., `NA` values). Valid options are:"
 )


}

roxy_na_action_fail <- function(data_label){
 paste0(
  "'fail' : an error is thrown if `", data_label,"` contains `NA` values"
 )
}

roxy_na_action_pass <- function(data_label){
 paste0(
  "'pass' : the output will have `NA` in all rows where `", data_label,"` has 1 or more `NA` value for the predictors used by `object`"
 )
}

roxy_na_action_omit <- function(data_label){
 paste0(
  "'omit' : rows in `", data_label,"` with incomplete data will be dropped"
 )
}

roxy_na_action_impute_meanmode <- function(data_label){
 paste0(
  "'impute_meanmode' : missing values for continuous and categorical variables in `", data_label,"` will be imputed using the mean and mode, respectively"
 )
}


# citations ---------------------------------------------------------------



roxy_cite <- function(authors,
                      title,
                      journal,
                      date,
                      number=NULL,
                      doi=NULL,
                      url=NULL){

 if(!is.null(number)) number <- paste0('; ', number)
 if(!is.null(doi)) doi <- paste(" DOI:", doi)
 if(!is.null(url)) url <- paste(" URL:", url)

 ending <- paste(c(doi, url), collapse = '.')

 paste0(
  authors, '. ',
  title, '. ',
  '*', journal, '* ',
  date,
  number, '.',
  ending
 )
}

roxy_cite_breiman_2001 <- function(){

 roxy_cite(
  authors = "Breiman L",
  title = "Random forests",
  journal = "Machine learning",
  date = "2001 Oct",
  number = "45(1):5-32",
  doi = "10.1023/A:1010933404324"
 )

}

roxy_cite_ishwaran_2008 <- function(){

 roxy_cite(
  authors = "Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS",
  title = "Random survival forests",
  journal = "Annals of applied statistics",
  date = "2008 Sep",
  number = "2(3):841-60",
  doi = "10.1214/08-AOAS169"
 )

}

roxy_cite_jaeger_2019 <- function(){

 roxy_cite(
  authors = "Jaeger BC, Long DL, Long DM, Sims M, Szychowski JM, Min YI, Mcclure LA, Howard G, Simon N",
  title = "Oblique random survival forests",
  journal = "Annals of applied statistics",
  date = "2019 Sep",
  number = "13(3):1847-83",
  doi = "10.1214/19-AOAS1261"
 )

}

roxy_cite_jaeger_2023 <- function(){

 roxy_cite(
  authors = "Jaeger BC, Welden S, Lenoir K, Speiser JL, Segar MW, Pandey A, Pajewski NM",
  title = "Accelerated and interpretable oblique random survival forests",
  journal = "Journal of Computational and Graphical Statistics",
  date = "Published online 08 Aug 2023",
  number = NULL,
  doi = "10.1080/10618600.2023.2231048"
  # url = "https://doi.org/10.1080/10618600.2023.2231048"
 )

}

roxy_cite_hooker_2021 <- function(){

 roxy_cite(
  authors = "Giles Hooker, Lucas Mentch, Siyu Zhou",
  title = "Unrestricted Permutation forces Extrapolation: Variable Importance Requires at least One More Model, or There Is No Free Variable Importance",
  journal = "arXiv e-prints",
  date = "2021 Oct",
  number = 'arXiv-1905',
  url = "https://doi.org/10.48550/arXiv.1905.03151"
 )

}

roxy_cite_harrell_1982 <- function(){

 roxy_cite(
  authors = "Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA",
  title = "Evaluating the Yield of Medical Tests",
  journal = "JAMA",
  date = "1982",
  number = '247(18):2543-2546',
  doi = "10.1001/jama.1982.03320430047030"
 )

}

roxy_cite_menze_2011 <- function(){

 roxy_cite(
  authors = "Menze BH, Kelm BM, Splitthoff DN, Koethe U, Hamprecht FA",
  title = "On oblique random forests",
  journal = "Joint European Conference on Machine Learning and Knowledge Discovery in Databases",
  date = "2011 Sep 4",
  number = 'pp. 453-469',
  doi = "10.1007/978-3-642-23783-6_29"
 )

}

roxy_cite_simon_2011 <- function(){

 roxy_cite(
  authors = "Simon N, Friedman J, Hastie T, Tibshirani R",
  title = "Regularization paths for Cox's proportional hazards model via coordinate descent",
  journal = "Journal of statistical software",
  date = "2011 Mar",
  number = '39(5):1',
  doi = "10.18637/jss.v039.i05"
 )

}


# ... ---------------------------------------------------------------------

roxy_dots <- function(){
 "Further arguments passed to or from other methods (not currently used)."
}



# variable importance descriptions ----------------------------------------

roxy_vi_describe <- function(type){

 switch(type,
        'negate' = "Each variable is assessed separately by multiplying the variable's coefficients by -1 and then determining how much the model's performance changes. The worse the model's performance after negating coefficients for a given variable, the more important the variable. This technique is promising b/c it does not require permutation and it emphasizes variables with larger coefficients in linear combinations, but it is also relatively new and hasn't been studied as much as permutation importance. See Jaeger, (2023) for more details on this technique.",
        'permute' = "Each variable is assessed separately by randomly permuting the variable's values and then determining how much the model's performance changes. The worse the model's performance after permuting the values of a given variable, the more important the variable. This technique is flexible, intuitive, and frequently used. It also has several [known limitations](https://christophm.github.io/interpretable-ml-book/feature-importance.html#disadvantages-9)",
        'anova' = "A p-value is computed for each coefficient in each linear combination of variables in each decision tree. Importance for an individual predictor variable is the proportion of times a p-value for its coefficient is < 0.01. This technique is very efficient computationally, but may not be as effective as permutation or negation in terms of selecting signal over noise variables. See [Menze, 2011](https://link.springer.com/chapter/10.1007/978-3-642-23783-6_29) for more details on this technique.")

}



# partial dependence ------------------------------------------------------

roxy_pd_oob_explain <- function(label){

 paste0("You can compute ", label, " ",
        "three ways using a random forest: \n",
        "- using in-bag predictions for the training data\n",
        "- using out-of-bag predictions for the training data\n",
        "- using predictions for a new set of data\n\n",
        "See examples for more details")

}

roxy_pd_explain <- function(){
 "Partial dependence (PD) shows the expected prediction from a model as a function of a single predictor or multiple predictors. The expectation is marginalized over the values of all other predictors, giving something like a multivariable adjusted estimate of the model's prediction."
}

roxy_ice_explain <- function(){
 "Unlike partial dependence, which shows the expected prediction as a function of one or multiple predictors, individual conditional expectations (ICE) show the prediction for an individual observation as a function of a predictor."
}

roxy_pd_limitations <- function(){
 "Partial dependence has a number of [known limitations and assumptions](https://christophm.github.io/interpretable-ml-book/pdp.html#disadvantages-5) that users should be aware of (see Hooker, 2021). In particular, partial dependence is less intuitive when >2 predictors are examined jointly, and it is assumed that the feature(s) for which the partial dependence is computed are not correlated with other features (this is likely not true in many cases). Accumulated local effect plots can be used (see [here](https://christophm.github.io/interpretable-ml-book/ale.html)) in the case where feature independence is not a valid assumption."
}

# nocov end
