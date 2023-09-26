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
                  RObject lincomb_R_function,
                  // predictions
                  PredType pred_type,
                  bool pred_mode,
                  bool pred_aggregate,
                  bool oobag_pred,
                  EvalType oobag_eval_type,
                  arma::uword oobag_eval_every,
                  Rcpp::RObject oobag_R_function,
                  uint n_thread,
                  int verbosity){

 this->data = std::move(input_data);
 this->tree_seeds = tree_seeds;
 this->n_tree = n_tree;
 this->mtry = mtry;
 this->grow_mode = grow_mode;
 this->vi_type = vi_type;
 this->vi_max_pvalue = vi_max_pvalue;
 this->leaf_min_obs = leaf_min_obs;
 this->split_rule = split_rule;
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
 this->pred_mode = pred_mode;
 this->pred_aggregate = pred_aggregate;
 this->oobag_pred = oobag_pred;
 this->oobag_eval_type = oobag_eval_type;
 this->oobag_eval_every = oobag_eval_every;
 this->oobag_R_function = oobag_R_function;
 this->n_thread = n_thread;
 this->verbosity = verbosity;

 if(vi_type != VI_NONE){
  vi_numer.zeros(data->get_n_cols());
  if(vi_type == VI_ANOVA){
   vi_denom.zeros(data->get_n_cols());
  }
 }

 // oobag denominator tracks the number of times an obs is oobag
 oobag_denom.zeros(data->get_n_rows());

 if(verbosity > 1){

  Rcout << "------------ input data dimensions ------------" << std::endl;
  Rcout << "N observations total: " << data->get_n_rows()    << std::endl;
  Rcout << "N columns total: "      << data->get_n_cols()    << std::endl;
  Rcout << "-----------------------------------------------";
  Rcout << std::endl;
  Rcout << std::endl;

 }

}

void Forest::run(bool oobag){

 if(pred_mode){

  init_trees();

  this->pred_values = predict(oobag);

 } else if (grow_mode) {

  // initialize the trees
  plant();
  // grow the trees
  grow();

  // compute out-of-bag predictions if needed
  if(oobag){
   this->pred_values = predict(oobag);
  }

 }

 if(vi_type == VI_PERMUTE || vi_type == VI_NEGATE){
  compute_oobag_vi();
 }

}

void Forest::init_trees(){

 for(uword i = 0; i < n_tree; ++i){

  trees[i]->init(data.get(),
                 tree_seeds[i],
                 mtry,
                 pred_type,
                 leaf_min_obs,
                 vi_type,
                 vi_max_pvalue,
                 split_rule,
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
                 lincomb_R_function,
                 verbosity);

 }

}

void Forest::grow() {

 // initialize trees before doing anything else
 init_trees();

 // Create thread ranges
 equalSplit(thread_ranges, 0, n_tree - 1, n_thread);

 // reset progress to 0
 progress = 0;

 if(n_thread == 1){
  // ensure safe usage of R functions and glmnet
  // by growing trees in a single thread.
  grow_single_thread(&vi_numer, &vi_denom);
  return;
 }

 // catch interrupts from threads
 aborted = false;
 aborted_threads = 0;

 // containers
 std::vector<std::thread> threads;
 std::vector<vec> vi_numer_threads(n_thread);
 std::vector<uvec> vi_denom_threads(n_thread);

 // reserve memory
 threads.reserve(n_thread);

 // begin multi-thread grow
 for (uint i = 0; i < n_thread; ++i) {

  vi_numer_threads[i].zeros(data->n_cols);
  if(vi_type == VI_ANOVA) vi_denom_threads[i].zeros(data->n_cols);

  threads.emplace_back(&Forest::grow_multi_thread, this, i,
                       &(vi_numer_threads[i]),
                       &(vi_denom_threads[i]));
 }

 if(verbosity == 1){
  show_progress("Growing trees", n_tree);
 }

 // end multi-thread grow
 for (auto &thread : threads) {
  thread.join();
 }

 if (aborted_threads > 0) {
  throw std::runtime_error("User interrupt.");
 }

 if(vi_type == VI_ANOVA){

  for(uint i = 0; i < n_thread; ++i){
   vi_numer += vi_numer_threads[i];
   vi_denom += vi_denom_threads[i];
  }

 }

}

void Forest::grow_single_thread(vec* vi_numer_ptr,
                                uvec* vi_denom_ptr){


 using std::chrono::steady_clock;
 using std::chrono::duration_cast;
 using std::chrono::seconds;

 steady_clock::time_point start_time = steady_clock::now();
 steady_clock::time_point last_time = steady_clock::now();
 size_t max_progress = n_tree;

 for (uint i = 0; i < n_tree; ++i) {

  if(verbosity > 1){
   Rcout << "------------ Growing tree " << i << " --------------";
   Rcout << std::endl;
   Rcout << std::endl;
  }

  trees[i]->grow(vi_numer_ptr, vi_denom_ptr);

  ++progress;

  if(verbosity == 1){

   seconds elapsed_time = duration_cast<seconds>(steady_clock::now() - last_time);

   if ((progress > 0 && elapsed_time.count() > STATUS_INTERVAL) ||
       (progress == max_progress)) {

    double relative_progress = (double) progress / (double) max_progress;
    seconds time_from_start = duration_cast<seconds>(steady_clock::now() - start_time);
    uint remaining_time = (1 / relative_progress - 1) * time_from_start.count();

    Rcout << "Growing trees: ";
    Rcout << round(100 * relative_progress) << "%. ";

    if(progress < max_progress){
     Rcout << "~ time remaining: ";
     Rcout << beautifyTime(remaining_time) << ".";
    }

    Rcout << std::endl;

    last_time = steady_clock::now();

   }

  }

  Rcpp::checkUserInterrupt();

 }

}


void Forest::grow_multi_thread(uint thread_idx,
                               vec* vi_numer_ptr,
                               uvec* vi_denom_ptr) {


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

void Forest::compute_oobag_vi() {

 // catch interrupts from threads
 aborted = false;
 aborted_threads = 0;

 // show progress from threads
 progress = 0;

 if(n_thread == 1){
  vec* vi_numer_ptr = &vi_numer;
  compute_oobag_vi_single_thread(vi_numer_ptr);
  return;
 }

 std::vector<std::thread> threads;
 std::vector<vec> vi_numer_threads(n_thread);
 // no denominator b/c it is equal to n_tree for all oob vi methods

 threads.reserve(n_thread);

 for (uint i = 0; i < n_thread; ++i) {

  vi_numer_threads[i].zeros(data->n_cols);

  threads.emplace_back(&Forest::compute_oobag_vi_multi_thread,
                       this, i, &(vi_numer_threads[i]));
 }

 if(verbosity == 1){
  show_progress("Computing importance", n_tree);
 }

 for (auto &thread : threads) {
  thread.join();
 }

 if (aborted_threads > 0) {
  throw std::runtime_error("User interrupt.");
 }

 for(uint i = 0; i < n_thread; ++i){
  vi_numer += vi_numer_threads[i];
 }

}

void Forest::compute_oobag_vi_single_thread(vec* vi_numer_ptr) {

 using std::chrono::steady_clock;
 using std::chrono::duration_cast;
 using std::chrono::seconds;

 steady_clock::time_point start_time = steady_clock::now();
 steady_clock::time_point last_time = steady_clock::now();
 size_t max_progress = n_tree;

 for(uint i = 0; i < n_tree; ++i){

  trees[i]->compute_oobag_vi(vi_numer_ptr, vi_type);

  ++progress;

  if(verbosity == 1){

   seconds elapsed_time = duration_cast<seconds>(steady_clock::now() - last_time);

   if ((progress > 0 && elapsed_time.count() > STATUS_INTERVAL) ||
       (progress == max_progress)) {

    double relative_progress = (double) progress / (double) max_progress;
    seconds time_from_start = duration_cast<seconds>(steady_clock::now() - start_time);
    uint remaining_time = (1 / relative_progress - 1) * time_from_start.count();

    Rcout << "Computing importance: ";
    Rcout << round(100 * relative_progress) << "%. ";

    if(progress < max_progress){
     Rcout << "~ time remaining: ";
     Rcout << beautifyTime(remaining_time) << ".";
    }

    Rcout << std::endl;

    last_time = steady_clock::now();

   }

  }

  Rcpp::checkUserInterrupt();

 }

}

void Forest::compute_oobag_vi_multi_thread(uint thread_idx, vec* vi_numer_ptr) {

 if (thread_ranges.size() > thread_idx + 1) {

  for(uint i=thread_ranges[thread_idx]; i<thread_ranges[thread_idx+1]; ++i){

   trees[i]->compute_oobag_vi(vi_numer_ptr, vi_type);

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

void Forest::compute_prediction_accuracy(Data* prediction_data,
                                         arma::mat& prediction_values,
                                         arma::uword row_fill){

 // avoid dividing by zero
 uvec valid_observations = find(oobag_denom > 0);

 // subset each data input
 mat y_valid = prediction_data->y_rows(valid_observations);
 vec w_valid = prediction_data->w_subvec(valid_observations);
 mat p_valid = prediction_values.rows(valid_observations);

 // scale predictions based on how many trees contributed
 // (important to note it's different for each oobag obs)
 vec valid_denom = oobag_denom(valid_observations);
 p_valid.each_col() /= valid_denom;

 // pass along to forest-specific version
 compute_prediction_accuracy(y_valid, w_valid, p_valid, row_fill);

}

mat Forest::predict(bool oobag) {

 mat result;

 // No. of cols in pred mat depend on the type of forest
 resize_pred_mat(result);

 // Slots to hold oobag prediction accuracy
 // (needs to be resized even if !oobag)
 resize_oobag_eval();

 progress = 0;
 aborted = false;
 aborted_threads = 0;

 if(n_thread == 1){
  // ensure safe usage of R functions
  predict_single_thread(data.get(), oobag, result);

 } else {

  std::vector<std::thread> threads;
  std::vector<mat> result_threads(n_thread);
  std::vector<vec> oobag_denom_threads(n_thread);

  threads.reserve(n_thread);

  for (uint i = 0; i < n_thread; ++i) {

   resize_pred_mat(result_threads[i]);
   if(oobag) oobag_denom_threads[i].zeros(data->n_rows);

   threads.emplace_back(&Forest::predict_multi_thread,
                        this, i, data.get(), oobag,
                        &(result_threads[i]),
                        &(oobag_denom_threads[i]));
  }

  if(verbosity == 1){
   show_progress("Computing predictions", n_tree);
  }

  // wait for all threads to finish before proceeding
  for (auto &thread : threads) {
   thread.join();
  }

  for(uint i = 0; i < n_thread; ++i){

   result += result_threads[i];

   if(oobag){

    oobag_denom += oobag_denom_threads[i];

    // evaluate oobag error after joining each thread
    // (only safe to do this when the condition below holds)
    if(n_tree/oobag_eval_every == n_thread && i<n_thread-1){

     // i should be uint to access threads,
     // eval_row should be uword to access oobag_eval
     uword eval_row = i;

     compute_prediction_accuracy(data.get(), result, eval_row);

    }
   }

  }

 }

 if(pred_type == PRED_TERMINAL_NODES || !pred_aggregate){
  return(result);
 }

 if(oobag){

  compute_prediction_accuracy(data.get(), result, oobag_eval.n_rows-1);

  result.each_col() /= oobag_denom;

 } else {

  result /= n_tree;

 }

 return(result);

}

void Forest::predict_single_thread(Data* prediction_data,
                                   bool oobag,
                                   mat& result) {

 using std::chrono::steady_clock;
 using std::chrono::duration_cast;
 using std::chrono::seconds;

 steady_clock::time_point start_time = steady_clock::now();
 steady_clock::time_point last_time = steady_clock::now();
 size_t max_progress = n_tree;

 for (uint i = 0; i < n_tree; ++i) {

  if(verbosity > 1){
   if(oobag){
    Rcout << "--- Computing oobag predictions: tree " << i << " ---";
   } else {
    Rcout << "------ Computing predictions: tree " << i << " -----";
   }
   Rcout << std::endl;
   Rcout << std::endl;
  }

  trees[i]->predict_leaf(prediction_data, oobag);

  if(pred_type == PRED_TERMINAL_NODES){
   result.col(i) = conv_to<vec>::from(trees[i]->get_pred_leaf());
  } else if (!pred_aggregate){
   vec col_i = result.unsafe_col(i);
   trees[i]->predict_value(&col_i, &oobag_denom, pred_type, oobag);
  } else {
   trees[i]->predict_value(&result, &oobag_denom, pred_type, oobag);
  }

  progress++;

  if(verbosity == 1){

   seconds elapsed_time = duration_cast<seconds>(steady_clock::now() - last_time);

   if ((progress > 0 && elapsed_time.count() > STATUS_INTERVAL) ||
       (progress == max_progress)) {

    double relative_progress = (double) progress / (double) max_progress;
    seconds time_from_start = duration_cast<seconds>(steady_clock::now() - start_time);
    uint remaining_time = (1 / relative_progress - 1) * time_from_start.count();

    Rcout << "Computing predictions: ";
    Rcout << round(100 * relative_progress) << "%. ";

    if(progress < max_progress){
     Rcout << "~ time remaining: ";
     Rcout << beautifyTime(remaining_time) << ".";
    }

    Rcout << std::endl;

    last_time = steady_clock::now();

   }

  }

  // if tracking oobag error over time:
  if(oobag && (progress % oobag_eval_every == 0) ){

   uword eval_row = (progress / oobag_eval_every) - 1;
   // mat preds = result.each_col() / oobag_denom;
   compute_prediction_accuracy(prediction_data, result, eval_row);

  }

 }

}

void Forest::predict_multi_thread(uint thread_idx,
                                  Data* prediction_data,
                                  bool oobag,
                                  mat* result_ptr,
                                  vec* denom_ptr) {

 if (thread_ranges.size() > thread_idx + 1) {

  for (uint i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {

   trees[i]->predict_leaf(prediction_data, oobag);

   trees[i]->predict_value(result_ptr, denom_ptr, pred_type, oobag);

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

arma::uword Forest::find_max_eval_steps(){

 if(!oobag_pred) return(0);

 uword n_evals = std::ceil(n_tree / oobag_eval_every);

 if(n_evals > n_tree) n_evals = n_tree;

 if(n_evals < 1) n_evals = 1;

 return(n_evals);

}

void Forest::resize_oobag_eval(){

 uword n_evals = find_max_eval_steps();

 oobag_eval.resize(n_evals, 1);

}

void Forest::show_progress(std::string operation, size_t max_progress) {

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

  if ((progress > 0 && elapsed_time.count() > STATUS_INTERVAL) ||
      (progress == max_progress)) {

   double relative_progress = (double) progress / (double) max_progress;
   seconds time_from_start = duration_cast<seconds>(steady_clock::now() - start_time);
   uint remaining_time = (1 / relative_progress - 1) * time_from_start.count();

   Rcout << operation << ": ";
   Rcout << round(100 * relative_progress) << "%. ";

   if(progress < max_progress){
    Rcout << "~ time remaining: ";
    Rcout << beautifyTime(remaining_time) << ".";
   }

   Rcout << std::endl;

   last_time = steady_clock::now();

  }
 }
}

}


