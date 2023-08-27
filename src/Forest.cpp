//  Forest.cpp

#include <RcppArmadillo.h>
#include "Forest.h"
#include "Tree.h"

using namespace arma;
using namespace Rcpp;

namespace aorsf {

Forest::Forest(){ }

void Forest::load(arma::uword n_tree,
                  std::vector<std::vector<double>>& forest_cutpoint,
                  std::vector<std::vector<arma::uword>>& forest_child_left,
                  std::vector<std::vector<arma::vec>>& forest_coef_values,
                  std::vector<std::vector<arma::uvec>>& forest_coef_indices,
                  std::vector<std::vector<arma::vec>>& forest_leaf_pred_horizon,
                  std::vector<std::vector<arma::vec>>& forest_leaf_pred_surv,
                  std::vector<std::vector<arma::vec>>& forest_leaf_pred_chf) {

 this->n_tree = n_tree;

 if(VERBOSITY > 0){
  Rcout << "---- loading forest from input list ----";
  Rcout << std::endl << std::endl;
 }


 // Create trees
 trees.reserve(n_tree);

 for (uword i = 0; i < n_tree; ++i) {
  trees.push_back(
   std::make_unique<Tree>(forest_cutpoint[i],
                          forest_child_left[i],
                          forest_coef_values[i],
                          forest_coef_indices[i],
                          forest_leaf_pred_horizon[i],
                          forest_leaf_pred_surv[i],
                          forest_leaf_pred_chf[i])
  );
 }

 // Create thread ranges
 equalSplit(thread_ranges, 0, n_tree - 1, n_thread);

}

void Forest::init(std::unique_ptr<Data> input_data,
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
                  uint n_thread){

 this->data = std::move(input_data);
 this->tree_seeds = tree_seeds;
 this->n_tree = n_tree;
 this->mtry = mtry;
 this->vi_type = vi_type;
 this->vi_max_pvalue = vi_max_pvalue;
 this->leaf_min_events = leaf_min_events;
 this->leaf_min_obs = leaf_min_obs;
 this->split_rule = split_rule;
 this->split_min_events = split_min_events;
 this->split_min_obs = split_min_obs;
 this->split_min_stat = split_min_stat;
 this->split_max_cuts = split_max_cuts;
 this->split_max_retry = split_max_retry;
 this->lincomb_type = lincomb_type; this->lincomb_eps = lincomb_eps;
 this->lincomb_iter_max = lincomb_iter_max;
 this->lincomb_scale = lincomb_scale;
 this->lincomb_alpha = lincomb_alpha;
 this->lincomb_df_target = lincomb_df_target;
 this->lincomb_ties_method = lincomb_ties_method;
 this->lincomb_R_function = lincomb_R_function;
 this->pred_type = pred_type;
 this->pred_horizon = pred_horizon;
 this->oobag_pred = oobag_pred;
 this->oobag_eval_every = oobag_eval_every;
 this->oobag_R_function = oobag_R_function;
 this->n_thread = n_thread;

 if(vi_type != VI_NONE){
  vi_numer.zeros(data->get_n_cols());
  vi_denom.zeros(data->get_n_cols());
 }

  if(VERBOSITY > 0){
  Rcout << "------------ input data dimensions ------------"   << std::endl;
  Rcout << "N obs total: "     << data->get_n_rows() << std::endl;
  Rcout << "N columns total: " << data->get_n_cols() << std::endl;
  Rcout << "-----------------------------------------------";
  Rcout << std::endl << std::endl;
 }

}

void Forest::init_trees(){

 for(uword i = 0; i < n_tree; ++i){

  trees[i]->init(data.get(),
                 tree_seeds[i],
                 mtry,
                 leaf_min_events,
                 leaf_min_obs,
                 vi_type,
                 vi_max_pvalue,
                 split_rule,
                 split_min_events,
                 split_min_obs,
                 split_min_stat,
                 split_max_cuts,
                 split_max_retry,
                 lincomb_type,
                 lincomb_eps,
                 lincomb_iter_max,
                 lincomb_scale,
                 lincomb_alpha,
                 lincomb_df_target,
                 lincomb_ties_method,
                 lincomb_R_function);

 }

}

// growInternal() in ranger
void Forest::plant() {

 trees.reserve(n_tree);

 for (arma::uword i = 0; i < n_tree; ++i) {
  trees.push_back(std::make_unique<Tree>());
 }

}

void Forest::grow() {

 init_trees();

 // if multi-threading isn't required
 if(n_thread == 1){

  vec* vi_numer_ptr = &this->vi_numer;
  uvec* vi_denom_ptr = &this->vi_denom;

  for(uint i = 0; i < n_tree; ++i){

   // Abort if user has pressed Ctrl + C or Escape in R.
   checkUserInterrupt();
   trees[i]->grow(vi_numer_ptr, vi_denom_ptr);

  }

  return;

 }

 // Create thread ranges
 equalSplit(thread_ranges, 0, n_tree - 1, n_thread);
 // catch interrupts from threads
 aborted = false;
 aborted_threads = 0;
 // show progress from threads
 progress = 0;

 std::vector<std::thread> threads;
 threads.reserve(n_thread);

 for (uint i = 0; i < n_thread; ++i) {
  threads.emplace_back(&Forest::grow_in_threads, this, i);
 }

 showProgress("Growing trees...", n_tree);

 for (auto &thread : threads) {
  thread.join();
 }

 if (aborted_threads > 0) {
  throw std::runtime_error("User interrupt.");
 }



}

void Forest::grow_in_threads(uint thread_idx) {

 vec* vi_numer_ptr = &this->vi_numer;
 uvec* vi_denom_ptr = &this->vi_denom;

 if (thread_ranges.size() > thread_idx + 1) {

  for (uint i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {

   trees[i]->grow(vi_numer_ptr, vi_denom_ptr);

   // Check for user interrupt
   if (aborted) {
    std::unique_lock<std::mutex> lock(mutex);
    ++aborted_threads;
    condition_variable.notify_one();
    return;
   }

   // Increase progress by 1 tree
   std::unique_lock<std::mutex> lock(mutex);
   ++progress;
   condition_variable.notify_one();

  }

 }

}

mat Forest::predict() {

 mat result(data->n_rows, pred_horizon.size());

 mat* result_ptr = &result;

 for(uint i = 0; i < n_tree; ++i){
  trees[i]->predict_leaf(data.get());
  trees[i]->predict_value(result_ptr, pred_horizon, 'S');
 }

 result /= n_tree;

 return(result);

}

void Forest::run() {

}

void Forest::showProgress(std::string operation, size_t max_progress) {

 using std::chrono::steady_clock;
 using std::chrono::duration_cast;
 using std::chrono::seconds;

 steady_clock::time_point start_time = steady_clock::now();
 steady_clock::time_point last_time = steady_clock::now();
 std::unique_lock<std::mutex> lock(mutex);

 // Wait for message from threads and show output if enough time elapsed
 while (progress < max_progress) {
  condition_variable.wait(lock);
  seconds elapsed_time = duration_cast<seconds>(steady_clock::now() - last_time);

  // Check for user interrupt
  if (!aborted && checkInterrupt()) {
   aborted = true;
  }
  if (aborted && aborted_threads >= n_thread) {
   return;
  }

  if (progress > 0 && elapsed_time.count() > STATUS_INTERVAL) {

   double relative_progress = (double) progress / (double) max_progress;
   seconds time_from_start = duration_cast<seconds>(steady_clock::now() - start_time);
   uint remaining_time = (1 / relative_progress - 1) * time_from_start.count();

   Rcout << operation << "Progress: ";
   Rcout << round(100 * relative_progress) << "%. ";
   Rcout << "Estimated remaining time: ";
   Rcout << beautifyTime(remaining_time) << ".";
   Rcout << std::endl;

   last_time = steady_clock::now();

  }
 }
}

}


