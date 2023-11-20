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
  arma::uword n_class,
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
 this->n_class   = n_class;
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
                                        n_class,
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

 double result = 0;

 if(oobag_eval_type == EVAL_R_FUNCTION){

  // initialize function from tree object
  // (Functions can't be stored in C++ classes, but Robjects can)
  Rcpp::Function f_oobag_eval = Rcpp::as<Rcpp::Function>(oobag_R_function);


  // go through all columns if multi-class y,
  // but only go through one column if y is binary
  // uword start = 0;
  // if(n_class == 2) start = 1;

  Rcpp::NumericVector w_ = Rcpp::wrap(w);

  for(uword i = 0; i < predictions.n_cols; ++i){

   vec y_i = y.unsafe_col(i);
   vec p_i = predictions.unsafe_col(i);

   Rcpp::NumericVector y_ = Rcpp::wrap(y_i);
   Rcpp::NumericVector p_ = Rcpp::wrap(p_i);

   Rcpp::NumericVector R_result = f_oobag_eval(y_, w_, p_);

   double result_addon = R_result[0];

   result += result_addon;

  }

  oobag_eval(row_fill, 0) = result / predictions.n_cols;

  return;

 }

 for(uword i = 0; i < predictions.n_cols; i++){
  vec y_i = y.unsafe_col(i);
  vec p_i = predictions.unsafe_col(i);
  result += compute_cstat_clsf(y_i, w, p_i);
 }

 oobag_eval(row_fill, 0) = result / predictions.n_cols;

}

// growInternal() in ranger
void ForestClassification::plant() {

 trees.reserve(n_tree);

 for (arma::uword i = 0; i < n_tree; ++i) {
  trees.push_back(std::make_unique<TreeClassification>(this->n_class));
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


