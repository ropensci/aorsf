/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "TreeClassification.h"
#include "Coxph.h"
#include "utility.h"
// #include "NodeSplitStats.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 TreeClassification::TreeClassification() { }

 TreeClassification::TreeClassification(arma::uword n_class){
  this->n_class = n_class;
  this->binary = n_class == 2;
 }

 TreeClassification::TreeClassification(arma::uword n_obs,
                                        arma::uword n_class,
                                        arma::uvec& rows_oobag,
                                        std::vector<double>& cutpoint,
                                        std::vector<arma::uword>& child_left,
                                        std::vector<arma::vec>& coef_values,
                                        std::vector<arma::uvec>& coef_indices,
                                        std::vector<arma::vec>& leaf_pred_prob,
                                        std::vector<double>& leaf_summary) :
 Tree(rows_oobag, cutpoint, child_left, coef_values, coef_indices, leaf_summary),
 leaf_pred_prob(leaf_pred_prob){

  this->n_class = n_class;
  this->binary = n_class == 2;
  find_rows_inbag(n_obs);

 }

 void TreeClassification::resize_leaves(arma::uword new_size){

  leaf_pred_prob.resize(new_size);
  leaf_summary.resize(new_size);

 }


 double TreeClassification::compute_split_score(){

  double result=0;

  switch (split_rule) {

  case SPLIT_GINI: {

   vec y_i = y_node.unsafe_col(y_col_split);
   result = compute_gini(y_i, w_node, g_node);

   // gini index: lower is better, so
   // transform to make consistent with other stats
   result = (result-1) * -1;

   break;
  }

  case SPLIT_CONCORD: {

   vec y_i = y_node.unsafe_col(y_col_split);
   result = compute_cstat_clsf(y_i, w_node, g_node);

   // if the split has good 'anti-prediction':
   if(result < 0.50){ result = 1 - result; }

   break;

  }

  default:
   stop("invalid split rule");
   break;

  }

  return(result);

 }

 void TreeClassification::sprout_leaf_internal(uword node_id){

  vec pred_prob = compute_pred_prob(y_node, w_node);
  leaf_pred_prob[node_id] = pred_prob;
  leaf_summary[node_id] = pred_prob.index_max();

 }

 arma::uword TreeClassification::predict_value_internal(
   arma::uvec& pred_leaf_sort,
   arma::mat& pred_output,
   PredType pred_type,
   bool oobag
 ){

  uword n_preds_made = 0;

  if(pred_type == PRED_PROBABILITY){

   for(auto& it : pred_leaf_sort){

    uword leaf_id = pred_leaf[it];

    // the stopping condition for oobag predictions
    if(leaf_id == max_nodes) break;

    // usual case: a prediction matrix with one column per class
    if(pred_output.n_cols > 1){
     pred_output.row(it) += leaf_pred_prob[leaf_id].t();
     n_preds_made++;
    }

    // unusual case: a single column matrix
    // this occurs when pred_aggregate is false in forest.cpp's predict
    if(pred_output.n_cols == 1){
     pred_output.at(it, 0) += leaf_pred_prob[leaf_id][1];
     n_preds_made++;
    }

   }

  } else if(pred_type == PRED_CLASS){

   for(auto& it : pred_leaf_sort){

    uword leaf_id = pred_leaf[it];
    if(leaf_id == max_nodes) break;

    // usual case: a prediction matrix with one column per class
    if(pred_output.n_cols > 1){
     pred_output.at(it, leaf_summary[leaf_id])++;
     n_preds_made++;
    }

    // unusual case: a single column matrix
    // this occurs when pred_aggregate is false in forest.cpp's predict
    if(pred_output.n_cols == 1){
     pred_output.at(it, 0) = leaf_summary[leaf_id];
     n_preds_made++;
    }

   }

  }

  return(n_preds_made);

 }

 arma::mat TreeClassification::glm_fit(){

  vec y_col = y_node.unsafe_col(y_col_split);

  mat out = logreg_fit(x_node,
                       y_col,
                       w_node,
                       lincomb_scale,
                       lincomb_eps,
                       lincomb_iter_max);

  return(out);

 }

 arma::mat TreeClassification::glmnet_fit(){

  arma::vec y_col = y_node.unsafe_col(y_col_split);

  NumericMatrix xx = wrap(x_node);
  NumericMatrix yy = wrap(y_col);
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

 arma::mat TreeClassification::user_fit(){

  vec y_col = y_node.unsafe_col(y_col_split);

  NumericMatrix xx = wrap(x_node);
  NumericMatrix yy = wrap(y_col);
  NumericVector ww = wrap(w_node);

  // initialize function from tree object
  // (Functions can't be stored in C++ classes, but RObjects can)
  Function f_beta = as<Function>(lincomb_R_function);

  NumericMatrix beta_R = f_beta(xx, yy, ww);

  mat beta = mat(beta_R.begin(), beta_R.nrow(), beta_R.ncol(), false);

  return(beta);

 }

 double TreeClassification::compute_prediction_accuracy_internal(
   arma::mat& preds
 ){

  double result = 0, denom = preds.n_cols;

  uword start = 0;

  if(binary){

   start = 1;
   denom = 1;

  }

  if (oobag_eval_type == EVAL_R_FUNCTION){

   // initialize function from tree object
   // (Functions can't be stored in C++ classes, but RObjects can)
   Function f_oobag_eval = as<Function>(oobag_R_function);

   NumericVector w_ = wrap(w_oobag);

   for(uword i = start; i < preds.n_cols; ++i){

    vec y_i = y_oobag.unsafe_col(i);
    vec p_i = preds.unsafe_col(i);

    NumericVector y_ = wrap(y_i);
    NumericVector p_ = wrap(p_i);
    NumericVector R_result = f_oobag_eval(y_, w_, p_);

    double result_addon = R_result[0];

    result += result_addon;

   }

   return(result / denom);

  }

  for(uword i = start; i < y_oobag.n_cols; i++){
   vec y_i = y_oobag.unsafe_col(i);
   vec p_i = preds.unsafe_col(i);
   result += compute_cstat_clsf(y_i, w_oobag, p_i);
  }

  return result / denom;

 }

 uword TreeClassification::find_safe_mtry(){

  if(binary){ return find_safe_mtry_binary(); }

  return find_safe_mtry_multiclass();


 }

 uword TreeClassification::find_safe_mtry_binary(){

  // conditions to split a column:
  //   >= 3 events per predictor
  //   >= 3 non-events per predictor

  double safer_mtry = mtry;
  double y_sum_ctrls = sum(y_node.col(0));
  double y_sum_cases = sum(y_node.col(1));

  if(verbosity > 3){
   Rcout << "   -- Y sums (unweighted): ";
   Rcout << y_sum_cases << " cases, ";
   Rcout << y_sum_ctrls << " controls" << std::endl;
  }

  splittable_y_cols.zeros(1);

  if(y_sum_cases >= 3 && y_sum_ctrls >= 3){

   splittable_y_cols[0] = 1;
   y_col_split = 1;

   double min_count = y_sum_cases;
   if(y_sum_cases > y_sum_ctrls) min_count = y_sum_ctrls;

   if(lincomb_type != LC_GLMNET){

    while (min_count / safer_mtry < 3){
     --safer_mtry;
    }

   }

   uword out = safer_mtry;

   return(out);

  }

  if(verbosity > 3){
   Rcout << "   -- No y columns are splittable";
   Rcout << std::endl << std::endl;
  }

  return 0;

 }
 uword TreeClassification::find_safe_mtry_multiclass(){

  // conditions to split a column:
  //   >= 3 events per predictor
  //   >= 3 non-events per predictor

  double safer_mtry = mtry;

  double n = y_node.n_rows;
  vec y_sum_cases = sum(y_node, 0).t();
  vec y_sum_ctrls = n - y_sum_cases;

  if(verbosity > 3){

   for(uword i = 0; i < y_sum_cases.size(); ++i){
    Rcout << "   -- For column " << i << ": ";
    Rcout << y_sum_cases[i] << " cases, ";
    Rcout << y_sum_ctrls[i] << " controls (unweighted)" << std::endl;
   }
  }

  splittable_y_cols.zeros(y_node.n_cols);

  uword counter = 0;

  for(uword i = 0; i < y_node.n_cols; ++i){

   if(y_sum_cases[i] >= 3 && y_sum_ctrls[i] >= 3){
    splittable_y_cols[counter] = i;
    counter++;
   }

  }

  splittable_y_cols.resize(counter);

  if(counter == 0){

   if(verbosity > 3){
    Rcout << "   -- No y columns are splittable" << std::endl << std::endl;
   }

   return counter;

  }

  if(verbosity > 3){
   for(auto &i : splittable_y_cols){
    Rcout << "   -- Y column " << i << " is splittable" << std::endl;
   }
  }

  uword best_count = 0;

  for(auto& ycol : splittable_y_cols){

   uword min_count;

   if(y_sum_cases[ycol] <= y_sum_ctrls[ycol]){
    min_count = y_sum_cases[ycol];
   } else {
    min_count = y_sum_ctrls[ycol];
   }

   if(min_count > best_count){
    y_col_split = ycol;
    best_count = min_count;
   }

  }

  if(verbosity > 3){
   Rcout << "   -- Most splittable Y column: " << y_col_split << std::endl;
  }

  // glmnet can handle higher dimension x,
  // but other methods probably cannot.
  if(lincomb_type != LC_GLMNET){

   while (best_count / safer_mtry < 3){
    --safer_mtry;
   }

  }

  uword out = safer_mtry;

  return(out);

 }

 uword TreeClassification::get_n_col_vi(){

  return(n_class);

 }

 void TreeClassification::predict_value_vi(mat& pred_values){

  for(uword i = 0; i < pred_values.n_rows; ++i){
   pred_values.row(i) = leaf_pred_prob[pred_leaf[i]].t();
  }

 }

 } // namespace aorsf

