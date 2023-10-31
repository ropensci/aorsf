//  Forest.cpp

#include <RcppArmadillo.h>
#include "ForestClassification.h"
#include "TreeClassification.h"

#include <memory>

using namespace arma;
using namespace Rcpp;

namespace aorsf {

ForestClassification::ForestClassification() { }

ForestClassification::ForestClassification(arma::uword n_class){

 this->n_class = n_class;

}

void ForestClassification::resize_pred_mat_internal(arma::mat& p){

 p.zeros(data->n_rows, this->n_class);

 if(verbosity > 3){
  Rcout << "   -- pred mat size: " << p.n_rows << " rows by ";
  Rcout << p.n_cols << " columns." << std::endl << std::endl;
 }

}

void ForestClassification::load(
  arma::uword n_tree,
  arma::uword n_obs,
  std::vector<arma::uvec>& forest_rows_oobag,
  std::vector<std::vector<double>>& forest_cutpoint,
  std::vector<std::vector<arma::uword>>& forest_child_left,
  std::vector<std::vector<arma::vec>>& forest_coef_values,
  std::vector<std::vector<arma::uvec>>& forest_coef_indices,
  std::vector<std::vector<arma::vec>>& forest_leaf_pred_prob,
  std::vector<std::vector<double>>& forest_leaf_summary,
  PartialDepType pd_type,
  std::vector<arma::mat>& pd_x_vals,
  std::vector<arma::uvec>& pd_x_cols,
  arma::vec& pd_probs
) {

 this->n_tree    = n_tree;
 this->pd_type   = pd_type;
 this->pd_x_vals = pd_x_vals;
 this->pd_x_cols = pd_x_cols;
 this->pd_probs  = pd_probs;

 if(verbosity > 2){
  Rcout << "---- loading forest from input list ----";
  Rcout << std::endl << std::endl;
 }

 // Create trees
 trees.reserve(n_tree);

 for (uword i = 0; i < n_tree; ++i) {
  trees.push_back(
   std::make_unique<TreeClassification>(n_obs,
                                        forest_rows_oobag[i],
                                        forest_cutpoint[i],
                                        forest_child_left[i],
                                        forest_coef_values[i],
                                        forest_coef_indices[i],
                                        forest_leaf_pred_prob[i],
                                        forest_leaf_summary[i])
  );
 }

 if(n_thread > 1){
  // Create thread ranges
  equalSplit(thread_ranges, 0, n_tree - 1, n_thread);
 }


}

void ForestClassification::compute_prediction_accuracy_internal(
  arma::mat& y,
  arma::vec& w,
  arma::mat& predictions,
  arma::uword row_fill
) {

 double cstat_sum = 0;

 vec y_0 = 1 - sum(y, 1);
 mat y_augment = join_horiz(y_0, y);

 for(uword i = 0; i < predictions.n_cols; i++){
  vec y_i = y_augment.unsafe_col(i);
  vec p_i = predictions.unsafe_col(i);
  cstat_sum += compute_cstat_clsf(y_i, w, p_i);
 }

 oobag_eval(row_fill, 0) = cstat_sum / predictions.n_cols;

}

// growInternal() in ranger
void ForestClassification::plant() {

 trees.reserve(n_tree);

 for (arma::uword i = 0; i < n_tree; ++i) {
  trees.push_back(std::make_unique<TreeClassification>());
 }

}

std::vector<std::vector<arma::vec>> ForestClassification::get_leaf_pred_prob() {

 std::vector<std::vector<arma::vec>> result;

 result.reserve(n_tree);

 for (auto& tree : trees) {
  auto& temp = dynamic_cast<TreeClassification&>(*tree);
  result.push_back(temp.get_leaf_pred_prob());
 }

 return result;

}

}


