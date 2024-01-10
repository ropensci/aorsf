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
       action, ". Default is 0, which allows a suitable",
       " number of threads to be used based on availability.",
       sep = "")
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
