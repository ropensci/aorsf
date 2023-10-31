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
    pred_output.row(it) = leaf_pred_prob[leaf_id].t();

    n_preds_made++;
    if(oobag) pred_denom[it]++;

   }

  } else if(pred_type == PRED_CLASS){

   for(auto& it : pred_leaf_sort){

    uword leaf_id = pred_leaf[it];
    if(leaf_id == max_nodes) break;
    pred_output.row(it) = leaf_summary[leaf_id];

    n_preds_made++;
    if(oobag) pred_denom[it]++;

   }

  }

  return(n_preds_made);

 }

 arma::mat TreeClassification::glm_fit(){

  std::uniform_int_distribution<uword> udist_ycol(0, y_node.n_cols - 1);
  uword j = udist_ycol(random_number_generator);
  vec y_col = y_node.unsafe_col(j);

  mat out = logreg_fit(x_node,
                       y_col,
                       w_node,
                       lincomb_scale,
                       lincomb_eps,
                       lincomb_iter_max);

  return(out);

 }

 double TreeClassification::compute_prediction_accuracy_internal(
   arma::vec& preds
 ){

  double cstat_sum = 0;

  for(uword i = 0; i < preds.n_cols; i++){
   vec y_i = y_oobag.unsafe_col(i);
   vec p_i = preds.unsafe_col(i);
   cstat_sum += compute_cstat_clsf(y_i, w_oobag, p_i);
  }

  return cstat_sum / preds.n_cols;

 }

 arma::uword TreeClassification::find_safe_mtry(){

  double safer_mtry = mtry;

  if(lincomb_type == LC_GLM){

   // Need 3:1 ratio of unweighted events:predictors

   double n = y_node.n_rows;
   vec y_sums = sum(y_node, 0).t();
   vec y_0 = { n - sum(y_sums) };
   y_sums = join_vert(y_0, y_sums);


   double n_events_max = max(y_sums);

   while(n_events_max / safer_mtry < 3){
    --safer_mtry;
    if(safer_mtry == 0) break;
   }

  }

  uword out = safer_mtry;

  return(out);

 }

 } // namespace aorsf

