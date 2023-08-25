
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
           // predictions
           PredType pred_type,
           bool pred_mode,
           double pred_horizon,
           bool oobag_pred,
           arma::uword oobag_eval_every,
           uint n_thread);

 // Grow or predict
 void run();

 void grow(Function lincomb_R_function);

 void grow_in_threads(uint thread_idx);

 void plant();

 void showProgress(std::string operation, size_t max_progress);


 std::vector<std::vector<arma::uvec>> get_coef_indices() {

  std::vector<std::vector<arma::uvec>> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_coef_indices());
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
 SplitRule split_rule;
 double split_min_events;
 double split_min_obs;
 double split_min_stat;
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

 // predictions
 PredType pred_type;
 double   pred_horizon;
 bool     oobag_pred;
 arma::uword oobag_eval_every;

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
