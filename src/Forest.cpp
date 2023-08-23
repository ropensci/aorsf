//  Forest.cpp

#include <RcppArmadillo.h>
#include "Forest.h"
#include "Tree.h"

using namespace arma;
using namespace Rcpp;

namespace aorsf {

Forest::Forest(){ }

void Forest::init(std::unique_ptr<Data> input_data,
                  Rcpp::IntegerVector& tree_seeds,
                  arma::uword n_tree,
                  arma::uword mtry,
                  VariableImportance vi_type,
                  double vi_max_pvalue,
                  double leaf_min_events,
                  double leaf_min_obs,
                  SplitRule split_rule,
                  double split_min_events,
                  double split_min_obs,
                  double split_min_stat,
                  arma::uword split_max_cuts,
                  arma::uword split_max_retry,
                  LinearCombo lincomb_type,
                  double lincomb_eps,
                  arma::uword lincomb_iter_max,
                  bool   lincomb_scale,
                  double lincomb_alpha,
                  arma::uword lincomb_df_target,
                  arma::uword lincomb_ties_method,
                  PredType pred_type,
                  bool   pred_mode,
                  double pred_horizon,
                  bool   oobag_pred,
                  arma::uword oobag_eval_every,
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
 this->pred_type = pred_type;
 this->pred_horizon = pred_horizon;
 this->oobag_pred = oobag_pred;
 this->oobag_eval_every = oobag_eval_every;
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

// growInternal() in ranger
void Forest::plant() {

 trees.reserve(n_tree);

 for (arma::uword i = 0; i < n_tree; ++i) {
  trees.push_back(std::make_unique<Tree>());
 }

}

void Forest::grow(Function& lincomb_R_function){

 // Create thread ranges
 equalSplit(thread_ranges, 0, n_tree - 1, n_thread);

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
                 lincomb_ties_method);

 }

 std::vector<std::thread> threads;
 threads.reserve(n_thread);

 // Initialize importance per thread
 std::vector<arma::vec> vi_numer_threads;
 std::vector<arma::uvec> vi_denom_threads;


 for (uint i = 0; i < n_thread; ++i) {

  // vi_numer_threads[i].zeros(data->get_n_cols());
  // vi_denom_threads[i].zeros(data->get_n_cols());

  threads.emplace_back(&Forest::grow_in_threads, this, i);

 }


 for (auto &thread : threads) {
  thread.join();
 }

}

void Forest::grow_in_threads(uint thread_idx) {

 // vec vi_numer(data->get_n_cols(), fill::zeros);
 // uvec vi_denom(data->get_n_cols(), fill::zeros);

 vec* vi_numer_ptr = &this->vi_numer;
 uvec* vi_denom_ptr = &this->vi_denom;

 if (thread_ranges.size() > thread_idx + 1) {

  for (uint i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {

   trees[i]->grow(vi_numer_ptr, vi_denom_ptr);

   // // Check for user interrupt
   // if (aborted) {
   //  std::unique_lock<std::mutex> lock(mutex);
   //  ++aborted_threads;
   //  condition_variable.notify_one();
   //  return;
   // }

   // Increase progress by 1 tree
   std::unique_lock<std::mutex> lock(mutex);
   // ++progress;
   condition_variable.notify_one();

  }

 }

 // if(VERBOSITY > 1){
 //
  // Rcout << "-- test VI numerator ---" << std::endl;
  // Rcout << vi_numer << std::endl << std::endl;
  // Rcout << "-- test VI denominator ---" << std::endl;
  // Rcout << vi_denom << std::endl << std::endl;
 //
 // }

}

void Forest::run(){ }

}

