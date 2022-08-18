
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

# oobag_fun ---------------------------------------------------------------

roxy_oobag_fun_header <- function(){

 "(_function_) to be used for evaluating out-of-bag prediction accuracy"

}

roxy_oobag_fun_default <- function(
  stat_label = "Harrell's C-statistic (1982)"
){
  paste("When `oobag_fun = NULL` (the default),",
        stat_label,
        "is used to evaluate accuracy.")
 }

roxy_oobag_fun_user <- function(){
 "if you use your own `oobag_fun` note the following:"
}

roxy_oobag_fun_inputs <- function(){
 "`oobag_fun` should have two inputs: `y_mat` and `s_vec`"
}

roxy_oobag_fun_ymat <- function(){
  "`y_mat` is a two column matrix with first column named 'time', second named 'status'"
}

roxy_oobag_fun_svec <- function(){
 "`s_vec` is a numeric vector containing predicted survival probabilities."
}


roxy_oobag_fun_return <- function(){
 "`oobag_fun` should return a numeric output of length 1"
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

roxy_cite_jaeger_2022 <- function(){

 roxy_cite(
  authors = "Jaeger BC, Welden S, Lenoir K, Speiser JL, Segar MW, Pandey A, Pajewski NM",
  title = "Accelerated and interpretable oblique random survival forests",
  journal = "arXiv e-prints",
  date = "2022 Aug",
  number = 'arXiv-2208',
  url = "https://arxiv.org/abs/2208.01129"
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
        'negate' = "Each variable is assessed separately by multiplying the variable's coefficients by -1 and then determining how much the model's performance changes. The worse the model's performance after negating coefficients for a given variable, the more important the variable.",
        'permute' = "Each variable is assessed separately by randomly permuting the variable's values and then determining how much the model's performance changes. The worse the model's performance after permuting the values of a given variable, the more important the variable.",
        'anova' = "A p-value is computed for each coefficient in each linear combination of variables in each decision tree. Importance for an individual predictor variable is the proportion of times a p-value for its coefficient is < 0.01.")

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

