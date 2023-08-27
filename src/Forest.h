
//  Forest.h

#ifndef Forest_H
#define Forest_H

#include "Data.h"
#include "globals.h"
#include "utility.h"
#include "Tree.h"

#include <thread>
#include <mutex>
#include <condition_variable>

namespace aorsf {

class Forest {

public:

 // Constructor

 Forest();

 // deleting the copy constructor
 Forest(const Forest&) = delete;
 // deleting the copy assignment operator
 Forest& operator=(const Forest&) = delete;

 // Methods

 void load(arma::uword n_tree,
           std::vector<std::vector<double>>& forest_cutpoint,
           std::vector<std::vector<arma::uword>>& forest_child_left,
           std::vector<std::vector<arma::vec>>& forest_coef_values,
           std::vector<std::vector<arma::uvec>>& forest_coef_indices,
           std::vector<std::vector<arma::vec>>& forest_leaf_pred_horizon,
           std::vector<std::vector<arma::vec>>& forest_leaf_pred_surv,
           std::vector<std::vector<arma::vec>>& forest_leaf_pred_chf);

 void init(std::unique_ptr<Data> input_data,
           Rcpp::IntegerVector& tree_seeds,
           arma::uword n_tree,
           arma::uword mtry,
           VariableImportance vi_type,
           double vi_max_pvalue,
           // leaves
           double leaf_min_events,
           double leaf_min_obs,
           // node splitting
           SplitRule split_rule,
           double split_min_events,
           double split_min_obs,
           double split_min_stat,
           arma::uword split_max_cuts,
           arma::uword split_max_retry,
           // linear combinations
           LinearCombo lincomb_type,
           double lincomb_eps,
           arma::uword lincomb_iter_max,
           bool lincomb_scale,
           double lincomb_alpha,
           arma::uword lincomb_df_target,
           arma::uword lincomb_ties_method,
           RObject lincomb_R_function,
           // predictions
           PredType pred_type,
           bool pred_mode,
           arma::vec pred_horizon,
           bool oobag_pred,
           arma::uword oobag_eval_every,
           Rcpp::RObject oobag_R_function,
           uint n_thread);

 // Grow or predict
 void run();

 void init_trees();

 void grow();

 void grow_in_threads(uint thread_idx);

 void plant();

 arma::mat predict();

 void showProgress(std::string operation, size_t max_progress);

 std::vector<std::vector<double>> get_cutpoint() {

  std::vector<std::vector<double>> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_cutpoint());
  }

  return result;

 }

 std::vector<std::vector<arma::uword>> get_child_left() {

  std::vector<std::vector<arma::uword>> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_child_left());
  }

  return result;

 }

 std::vector<std::vector<arma::vec>> get_coef_values() {

  std::vector<std::vector<arma::vec>> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_coef_values());
  }

  return result;

 }
 std::vector<std::vector<arma::uvec>> get_coef_indices() {

  std::vector<std::vector<arma::uvec>> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_coef_indices());
  }

  return result;

 }
 std::vector<std::vector<arma::vec>> get_leaf_pred_horizon() {

  std::vector<std::vector<arma::vec>> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_leaf_pred_horizon());
  }

  return result;

 }
 std::vector<std::vector<arma::vec>> get_leaf_pred_surv() {

  std::vector<std::vector<arma::vec>> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_leaf_pred_surv());
  }

  return result;

 }
 std::vector<std::vector<arma::vec>> get_leaf_pred_chf() {

  std::vector<std::vector<arma::vec>> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_leaf_pred_chf());
  }

  return result;

 }


 // Member variables

 arma::uword n_tree;
 arma::uword mtry;
 Rcpp::IntegerVector tree_seeds;

 std::vector<std::unique_ptr<Tree>> trees;

 std::unique_ptr<Data> data;

 // variable importance
 VariableImportance vi_type;
 double vi_max_pvalue;

 arma::vec vi_numer;
 arma::uvec vi_denom;

 // leaves
 double leaf_min_events;
 double leaf_min_obs;

 // node splitting
 SplitRule   split_rule;
 double      split_min_events;
 double      split_min_obs;
 double      split_min_stat;
 arma::uword split_max_cuts;
 arma::uword split_max_retry;

 // linear combinations
 LinearCombo lincomb_type;
 double      lincomb_eps;
 bool        lincomb_scale;
 double      lincomb_alpha;
 arma::uword lincomb_iter_max;
 arma::uword lincomb_df_target;
 arma::uword lincomb_ties_method;
 RObject     lincomb_R_function;

 // predictions
 PredType pred_type;
 arma::vec pred_horizon;

 // out-of-bag
 bool        oobag_pred;
 arma::uword oobag_eval_every;
 RObject     oobag_R_function;


 // multi-threading
 uint n_thread;
 std::vector<uint> thread_ranges;
 std::mutex mutex;
 std::condition_variable condition_variable;

 size_t progress;
 size_t aborted_threads;
 bool aborted;


};

}



#endif /* Forest_H */
