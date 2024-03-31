/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "TreeRegression.h"
#include "Coxph.h"
#include "utility.h"
// #include "NodeSplitStats.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 TreeRegression::TreeRegression() { }

 TreeRegression::TreeRegression(arma::uword n_obs,
                                arma::uvec& rows_oobag,
                                std::vector<double>& cutpoint,
                                std::vector<arma::uword>& child_left,
                                std::vector<arma::vec>& coef_values,
                                std::vector<arma::uvec>& coef_indices,
                                std::vector<arma::vec>& leaf_pred_prob,
                                std::vector<double>& leaf_summary) :
 Tree(rows_oobag, cutpoint, child_left, coef_values, coef_indices, leaf_summary),
 leaf_pred_prob(leaf_pred_prob){

  find_rows_inbag(n_obs);

 }

 void TreeRegression::resize_leaves(arma::uword new_size){

  leaf_pred_prob.resize(new_size);
  leaf_summary.resize(new_size);

 }


 double TreeRegression::compute_split_score(){

  double result=0;

  switch (split_rule) {

  case SPLIT_VARIANCE: {

   for(uword i = 0; i < y_node.n_cols; i++){
   vec y_i = y_node.unsafe_col(i);
   result += compute_var_reduction(y_i, w_node, g_node);
  }

   result /= y_node.n_cols;
   break;

  }

  default:
   Rcpp::stop("invalid split rule");
   break;

  }

  return(result);

 }

 bool TreeRegression::is_node_splittable_internal(){

  // if y_node has < 3 unique values, you are done!

  double y_reference = y_node.at(0, 0);
  uword n_unique = 1;

  for(uword i = 0; i < y_node.n_rows; ++i){

   if(y_node.at(i, 0) != y_reference){
    n_unique++;
    y_reference = y_node.at(i, 0);
   }

   if(n_unique == 3) break;

  }

  double n_obs = sum(w_node);

  return(n_obs >= 2*leaf_min_obs &&
         n_obs >= split_min_obs &&
         n_unique >= 3);

 }

 void TreeRegression::sprout_leaf_internal(uword node_id){

  double pred_mean = compute_pred_mean(y_node, w_node);

  leaf_summary[node_id] = pred_mean;

  // vec quant_probs = {0.25, 0.50, 0.75};

  // // TODO: make weighted version
  // vec quant_vals = quantile(y_node, quant_probs);

  // leaf_pred_prob[node_id] = quant_vals;

 }

 arma::uword TreeRegression::predict_value_internal(
   arma::uvec& pred_leaf_sort,
   arma::mat& pred_output,
   PredType pred_type,
   bool oobag
 ){

  uword n_preds_made = 0;

  if(pred_type == PRED_PROBABILITY){

   for(auto& it : pred_leaf_sort){

    uword leaf_id = pred_leaf[it];
    if(leaf_id == max_nodes) break;
    pred_output.row(it) += leaf_pred_prob[leaf_id].t();

    n_preds_made++;

   }

  } else if(pred_type == PRED_MEAN){

   for(auto& it : pred_leaf_sort){

    uword leaf_id = pred_leaf[it];
    if(leaf_id == max_nodes) break;

    pred_output.at(it, 0) += leaf_summary[leaf_id];

    n_preds_made++;

   }

  }

  return(n_preds_made);

 }

 arma::mat TreeRegression::glm_fit(){

  vec y_col = y_node.unsafe_col(0);

  mat out = linreg_fit(x_node,
                       y_col,
                       w_node,
                       lincomb_scale,
                       lincomb_eps,
                       lincomb_iter_max);

  return(out);

 }

 arma::mat TreeRegression::glmnet_fit(){

  NumericMatrix xx = wrap(x_node);
  NumericMatrix yy = wrap(y_node);
  NumericVector ww = wrap(w_node);

  // initialize function from tree object
  // (Functions can't be stored in C++ classes, but RObjects can)
  Function f_beta = as<Function>(lincomb_R_function);

  NumericMatrix beta_R = f_beta(xx, yy, ww,
                                lincomb_alpha,
                                lincomb_df_target);

  mat beta = mat(beta_R.begin(), beta_R.nrow(), beta_R.ncol(), false);

  return(beta);

 }

 arma::mat TreeRegression::user_fit(){

  NumericMatrix xx = wrap(x_node);
  NumericMatrix yy = wrap(y_node);
  NumericVector ww = wrap(w_node);

  // initialize function from tree object
  // (Functions can't be stored in C++ classes, but RObjects can)
  Function f_beta = as<Function>(lincomb_R_function);

  NumericMatrix beta_R = f_beta(xx, yy, ww);

  mat beta = mat(beta_R.begin(), beta_R.nrow(), beta_R.ncol(), false);

  return(beta);

 }


 double TreeRegression::compute_prediction_accuracy_internal(
   arma::mat& preds
 ){

  if (oobag_eval_type == EVAL_R_FUNCTION){

   vec preds_vec = preds.unsafe_col(0);

   NumericMatrix y_wrap = wrap(y_oobag);
   NumericVector w_wrap = wrap(w_oobag);
   NumericVector p_wrap = wrap(preds_vec);

   // initialize function from tree object
   // (Functions can't be stored in C++ classes, but RObjects can)
   Function f_oobag = as<Function>(oobag_R_function);

   NumericVector result_R = f_oobag(y_wrap, w_wrap, p_wrap);

   return(result_R[0]);

  }

  double result = 0;

  for(uword i = 0; i < y_oobag.n_cols; i++){
   vec y_i = y_oobag.unsafe_col(i);
   vec p_i = preds.unsafe_col(i);
   result += compute_rsq(y_i, w_oobag, p_i);
  }

  return result / preds.n_cols;

 }

 arma::uword TreeRegression::find_safe_mtry(){

  double safer_mtry = mtry;

  if(lincomb_type == LC_GLM ||
     lincomb_type == LC_GLMNET){

   // conditions to split a column:
   //   >= 3 non-events per predictor

   double n = y_node.n_rows;

   if(verbosity > 3){
    Rcout << "   -- N obs (unweighted): " << n << std::endl;
   }

   while (n / safer_mtry < 3){
    --safer_mtry;
    if(safer_mtry == 0) break;
   }

  }

  uword out = safer_mtry;

  return(out);

 }

 uword TreeRegression::get_n_col_vi(){
  return(1);
 }

 void TreeRegression::predict_value_vi(mat& pred_values){

  for(uword i = 0; i < pred_values.n_rows; ++i){
   pred_values.at(i, 0) = leaf_summary[pred_leaf[i]];
  }

 }

 } // namespace aorsf

