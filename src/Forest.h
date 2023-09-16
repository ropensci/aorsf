
//  Forest.h

#ifndef FOREST_H
#define FOREST_H

#include "Data.h"
#include "globals.h"
#include "utility.h"
#include "Tree.h"
#include "TreeSurvival.h"

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

 virtual ~Forest() = default;

 // Methods

 void init(std::unique_ptr<Data> input_data,
           Rcpp::IntegerVector& tree_seeds,
           arma::uword n_tree,
           arma::uword mtry,
           VariableImportance vi_type,
           double vi_max_pvalue,
           // leaves
           double leaf_min_obs,
           // node splitting
           SplitRule split_rule,
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
           bool oobag_pred,
           arma::uword oobag_eval_every,
           Rcpp::RObject oobag_R_function,
           uint n_thread);

 // Grow or predict
 void run();

 std::vector<std::vector<double>> get_cutpoint() {

  std::vector<std::vector<double>> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_cutpoint());
  }

  return result;

 }
 std::vector<arma::uvec> get_rows_oobag() {

  std::vector<arma::uvec> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_rows_oobag());
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

 std::vector<std::vector<double>> get_leaf_summary() {

  std::vector<std::vector<double>> result;

  result.reserve(n_tree);

  for (auto& tree : trees) {
   result.push_back(tree->get_leaf_summary());
  }

  return result;

 }

 void set_unique_event_times(arma::vec& x){
  this->unique_event_times = x;
 }

 arma::vec& get_unique_event_times(){
  return(unique_event_times);
 }

 arma::vec& get_vi_numer(){
  return(vi_numer);
 }

 arma::uvec& get_vi_denom(){
  return(vi_denom);
 }

 virtual void plant() = 0;

 void grow();

 arma::mat predict(bool oobag);

protected:

 void init_trees();


 void grow_in_threads(uint thread_idx,
                      vec* vi_numer_ptr,
                      uvec* vi_denom_ptr);



 void predict_in_threads(uint thread_idx,
                         Data* prediction_data,
                         bool oobag,
                         mat* result_ptr,
                         vec* denom_ptr);

 void showProgress(std::string operation, size_t max_progress);

 virtual void resize_pred_mat(arma::mat& p) = 0;

 // Member variables

 arma::uword n_tree;
 arma::uword mtry;
 Rcpp::IntegerVector tree_seeds;

 std::vector<std::unique_ptr<Tree>> trees;

 std::unique_ptr<Data> data;

 arma::vec unique_event_times;

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
