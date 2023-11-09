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
 }

 TreeClassification::TreeClassification(arma::uword n_obs,
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

 void TreeClassification::resize_leaves(arma::uword new_size){

  leaf_pred_prob.resize(new_size);
  leaf_summary.resize(new_size);

 }


 double TreeClassification::compute_split_score(){

  double result=0;

  switch (split_rule) {

  case SPLIT_GINI: {

   for(uword i = 0; i < y_node.n_cols; i++){
   vec y_i = y_node.unsafe_col(i);
   result += compute_gini(y_i, w_node, g_node);
  }
   result /= y_node.n_cols;
   // gini index: lower is better, so
   // transform to make consistent with other stats
   result = (result-1) * -1;

   break;
  }

  case SPLIT_CONCORD: {

   for(uword i = 0; i < y_node.n_cols; i++){
    vec y_i = y_node.unsafe_col(i);
    result += compute_cstat_clsf(y_i, w_node, g_node);
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

 void TreeClassification::sprout_leaf_internal(uword node_id){

  vec pred_prob = compute_pred_prob(y_node, w_node);
  leaf_pred_prob[node_id] = pred_prob;
  leaf_summary[node_id] = pred_prob.index_max();

 }

 arma::uword TreeClassification::predict_value_internal(
   arma::uvec& pred_leaf_sort,
   arma::mat& pred_output,
   arma::vec& pred_denom,
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
    if(oobag) pred_denom[it]++;

   }

  } else if(pred_type == PRED_CLASS){

   for(auto& it : pred_leaf_sort){

    uword leaf_id = pred_leaf[it];
    if(leaf_id == max_nodes) break;

    pred_output.at(it, leaf_summary[leaf_id])++;

    // pred_output.row(it) = leaf_summary[leaf_id];

    n_preds_made++;
    if(oobag) pred_denom[it]++;

   }

  }

  return(n_preds_made);

 }

 arma::mat TreeClassification::glm_fit(){

  vec y_col;

  if(splittable_y_cols.size() > 1){
   std::uniform_int_distribution<uword> udist_ycol(0, splittable_y_cols.size() - 1);
   uword j = udist_ycol(random_number_generator);
   uword k = splittable_y_cols[j];
   y_col   = y_node.unsafe_col(k);
  } else {
   y_col   = y_node.unsafe_col(0);
  }


  mat out = logreg_fit(x_node,
                       y_col,
                       w_node,
                       lincomb_scale,
                       lincomb_eps,
                       lincomb_iter_max);

  return(out);

 }

 double TreeClassification::compute_prediction_accuracy_internal(
   arma::mat& preds
 ){

  double cstat_sum = 0;

  // note: preds includes a column for the non-case, but y does not.
  // That is why the preds column is ahead by 1 here.

  for(uword i = 0; i < y_oobag.n_cols; i++){
   vec y_i = y_oobag.unsafe_col(i);
   vec p_i = preds.unsafe_col(i+1);
   cstat_sum += compute_cstat_clsf(y_i, w_oobag, p_i);
  }

  return cstat_sum / preds.n_cols;

 }

 arma::uword TreeClassification::find_safe_mtry(){

  double safer_mtry = mtry;

  if(lincomb_type == LC_GLM){

   // conditions to split a column:
   //   >= 3 events per predictor
   //   >= 3 non-events per predictor

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

   for (auto& i : splittable_y_cols){

    while (y_sum_cases[i] / safer_mtry < 3 ||
           y_sum_ctrls[i] / safer_mtry < 3){
     --safer_mtry;
    }

   }

  }

  uword out = safer_mtry;

  return(out);

 }

 uword TreeClassification::get_n_col_vi(){
  return(n_class);
 }

 void TreeClassification::fill_pred_values_vi(mat& pred_values){

  for(uword i = 0; i < pred_values.n_rows; ++i){
   pred_values.row(i) = leaf_pred_prob[pred_leaf[i]].t();
  }

 }

 } // namespace aorsf

