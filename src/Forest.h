
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
           bool sample_with_replacement,
           double sample_fraction,
           bool grow_mode,
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
           Rcpp::RObject lincomb_R_function,
           // predictions
           PredType pred_type,
           bool pred_mode,
           bool pred_aggregate,
           PartialDepType pd_type,
           std::vector<arma::mat>& pd_x_vals,
           std::vector<arma::uvec>& pd_x_cols,
           arma::vec& pd_probs,
           bool oobag_pred,
           EvalType oobag_eval_type,
           arma::uword oobag_eval_every,
           Rcpp::RObject oobag_R_function,
           uint n_thread,
           int verbosity);

 // Grow or predict
 // void run(bool verbose, bool oobag);

 virtual void compute_prediction_accuracy(
   Data*       prediction_data,
   arma::mat&  prediction_values,
   arma::uword row_fill
 );

 virtual void compute_prediction_accuracy(
   arma::mat& y,
   arma::vec& w,
   arma::mat& predictions,
   arma::uword row_fill
 ) = 0;

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

 arma::mat& get_oobag_eval(){
  return(oobag_eval);
 }

 arma::mat& get_predictions(){
  return(pred_values);
 }

 std::vector<std::vector<arma::mat>>& get_pd_values(){
  return(pd_values);
 }

 void run(bool oobag);

 virtual void plant() = 0;

 void grow();

 arma::mat predict(bool oobag);

 std::vector<std::vector<arma::mat>> compute_dependence(bool oobag);

protected:

 void init_trees();

 void grow_single_thread(vec* vi_numer_ptr,
                         uvec* vi_denom_ptr);

 void grow_multi_thread(uint thread_idx,
                        vec* vi_numer_ptr,
                        uvec* vi_denom_ptr);

 void predict_single_thread(Data* prediction_data,
                            bool oobag,
                            mat& result);

 void predict_multi_thread(uint thread_idx,
                           Data* prediction_data,
                           bool oobag,
                           mat* result_ptr,
                           vec* denom_ptr);

 void compute_oobag_vi();

 void compute_oobag_vi_single_thread(vec* vi_numer_ptr);

 void compute_oobag_vi_multi_thread(uint thread_idx, vec* vi_numer_ptr);

 void show_progress(std::string operation, size_t max_progress);

 virtual void resize_pred_mat(arma::mat& p);

 virtual void resize_pred_mat_internal(arma::mat& p) = 0;

 arma::uword find_max_eval_steps();

 virtual void resize_oobag_eval();

 // Member variables

 arma::uword n_tree;
 arma::uword mtry;
 bool sample_with_replacement;
 double sample_fraction;
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
 Rcpp::RObject     lincomb_R_function;

 bool grow_mode;

 // predictions
 PredType pred_type;
 bool pred_mode;
 bool pred_aggregate;
 arma::mat pred_values;

 // partial dependence
 PartialDepType pd_type;
 std::vector<std::vector<arma::mat>> pd_values;
 std::vector<arma::mat> pd_x_vals;
 std::vector<arma::uvec> pd_x_cols;
 arma::vec pd_probs;

 // out-of-bag
 bool        oobag_pred;
 arma::vec   oobag_denom;
 arma::mat   oobag_eval;
 EvalType    oobag_eval_type;
 arma::uword oobag_eval_every;
 Rcpp::RObject     oobag_R_function;


 // multi-threading
 uint n_thread;
 std::vector<uint> thread_ranges;
 std::mutex mutex;
 std::condition_variable condition_variable;

 size_t progress;
 size_t aborted_threads;
 bool aborted;

 // printing to console
 int verbosity;


};

}



#endif /* Forest_H */
