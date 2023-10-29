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
   result = compute_gini(y_node, w_node, g_node);
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

 } // namespace aorsf

