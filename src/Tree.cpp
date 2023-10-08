/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "Tree.h"
#include "Coxph.h"

#include <memory>
#include <random>

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 Tree::Tree() :
   data(0),
   n_cols_total(0),
   n_rows_total(0),
   seed(0),
   mtry(0),
   pred_type(DEFAULT_PRED_TYPE),
   vi_type(VI_NONE),
   vi_max_pvalue(DEFAULT_ANOVA_VI_PVALUE),
   // leaf_min_events(DEFAULT_LEAF_MIN_EVENTS),
   leaf_min_obs(DEFAULT_LEAF_MIN_OBS),
   split_rule(DEFAULT_SPLITRULE),
   // split_min_events(DEFAULT_SPLIT_MIN_EVENTS),
   split_min_obs(DEFAULT_SPLIT_MIN_OBS),
   split_min_stat(DEFAULT_SPLIT_MIN_STAT),
   split_max_cuts(DEFAULT_SPLIT_MAX_CUTS),
   split_max_retry(DEFAULT_SPLIT_MAX_RETRY),
   lincomb_type(DEFAULT_LINCOMB),
   lincomb_eps(DEFAULT_LINCOMB_EPS),
   lincomb_iter_max(DEFAULT_LINCOMB_ITER_MAX),
   lincomb_scale(DEFAULT_LINCOMB_SCALE),
   lincomb_alpha(DEFAULT_LINCOMB_ALPHA),
   lincomb_df_target(0),
   lincomb_ties_method(DEFAULT_LINCOMB_TIES_METHOD),
   lincomb_R_function(0),
   verbosity(0){

 }

 Tree::Tree(arma::uvec& rows_oobag,
            std::vector<double>& cutpoint,
            std::vector<arma::uword>& child_left,
            std::vector<arma::vec>& coef_values,
            std::vector<arma::uvec>& coef_indices,
            std::vector<double>& leaf_summary) :
 data(0),
 n_cols_total(0),
 n_rows_total(0),
 seed(0),
 mtry(0),
 pred_type(DEFAULT_PRED_TYPE),
 vi_type(VI_NONE),
 vi_max_pvalue(DEFAULT_ANOVA_VI_PVALUE),
 // leaf_min_events(DEFAULT_LEAF_MIN_EVENTS),
 leaf_min_obs(DEFAULT_LEAF_MIN_OBS),
 split_rule(DEFAULT_SPLITRULE),
 // split_min_events(DEFAULT_SPLIT_MIN_EVENTS),
 split_min_obs(DEFAULT_SPLIT_MIN_OBS),
 split_min_stat(DEFAULT_SPLIT_MIN_STAT),
 split_max_cuts(DEFAULT_SPLIT_MAX_CUTS),
 split_max_retry(DEFAULT_SPLIT_MAX_RETRY),
 lincomb_type(DEFAULT_LINCOMB),
 lincomb_eps(DEFAULT_LINCOMB_EPS),
 lincomb_iter_max(DEFAULT_LINCOMB_ITER_MAX),
 lincomb_scale(DEFAULT_LINCOMB_SCALE),
 lincomb_alpha(DEFAULT_LINCOMB_ALPHA),
 lincomb_df_target(0),
 lincomb_ties_method(DEFAULT_LINCOMB_TIES_METHOD),
 lincomb_R_function(0),
 verbosity(0),
 rows_oobag(rows_oobag),
 cutpoint(cutpoint),
 child_left(child_left),
 coef_values(coef_values),
 coef_indices(coef_indices),
 leaf_summary(leaf_summary){

  this->max_nodes = cutpoint.size()+1;
  this->max_leaves = cutpoint.size()+1;

 }


 void Tree::init(Data* data,
                 int seed,
                 arma::uword mtry,
                 bool sample_with_replacement,
                 double sample_fraction,
                 PredType pred_type,
                 // double leaf_min_events,
                 double leaf_min_obs,
                 VariableImportance vi_type,
                 double vi_max_pvalue,
                 SplitRule split_rule,
                 // double split_min_events,
                 double split_min_obs,
                 double split_min_stat,
                 arma::uword split_max_cuts,
                 arma::uword split_max_retry,
                 LinearCombo lincomb_type,
                 double lincomb_eps,
                 arma::uword lincomb_iter_max,
                 bool lincomb_scale,
                 double lincomb_alpha,
                 arma::uword lincomb_df_target,
                 arma::uword lincomb_ties_method,
                 RObject lincomb_R_function,
                 RObject oobag_R_function,
                 EvalType oobag_eval_type,
                 int verbosity){

  // Initialize random number generator and set seed
  random_number_generator.seed(seed);

  this->data = data;
  this->n_cols_total = data->n_cols;
  this->n_rows_total = data->n_rows;
  this->seed = seed;
  this->mtry = mtry;
  this->sample_with_replacement = sample_with_replacement;
  this->sample_fraction = sample_fraction;
  this->pred_type = pred_type;
  // this->leaf_min_events = leaf_min_events;
  this->leaf_min_obs = leaf_min_obs;
  this->vi_type = vi_type;
  this->vi_max_pvalue = vi_max_pvalue;
  this->split_rule = split_rule;
  // this->split_min_events = split_min_events;
  this->split_min_obs = split_min_obs;
  this->split_min_stat = split_min_stat;
  this->split_max_cuts = split_max_cuts;
  this->split_max_retry = split_max_retry;
  this->lincomb_type = lincomb_type;
  this->lincomb_eps = lincomb_eps;
  this->lincomb_iter_max = lincomb_iter_max;
  this->lincomb_scale = lincomb_scale;
  this->lincomb_alpha = lincomb_alpha;
  this->lincomb_df_target = lincomb_df_target;
  this->lincomb_ties_method = lincomb_ties_method;
  this->lincomb_R_function = lincomb_R_function;
  this->oobag_R_function = oobag_R_function;
  this->oobag_eval_type = oobag_eval_type;
  this->verbosity = verbosity;

 }

 void Tree::resize_oobag_data(){

  if(rows_oobag.size() == 0){
   stop("attempting to allocate oob memory with empty rows_oobag");
  }

  x_oobag = data->x_rows(rows_oobag);
  y_oobag = data->y_rows(rows_oobag);
  w_oobag = data->w_subvec(rows_oobag);

 }

 // not currently used but will be in the future
 // # nocov start
 void Tree::resize_leaves(arma::uword new_size){

  leaf_summary.resize(new_size);

 }
 // # nocov end

 void Tree::sample_rows(){

  uword i, draw, n = data->n_rows;

  // Start with all samples OOB
  vec w_inbag(n, fill::zeros);

  std::uniform_int_distribution<uword> udist_rows(0, n - 1);

  if(sample_with_replacement){

   for (i = 0; i < n; ++i) {
    draw = udist_rows(random_number_generator);
    ++w_inbag[draw];
   }

  } else {

   if(sample_fraction == 1){
    w_inbag.fill(1);
   } else {

    uword n_sample = (uword) std::round(n * sample_fraction);
    for (i = 0; i < n_sample; ++i) {
     draw = udist_rows(random_number_generator);
     while(w_inbag[draw] == 1){
      draw = udist_rows(random_number_generator);
     }
     ++w_inbag[draw];

    }
   }

  }

  // multiply w_inbag by user specified weights.
  if(data->has_weights){
   w_inbag = w_inbag % data->w;
  }

  this->rows_inbag = find(w_inbag > 0);
  this->rows_oobag = find(w_inbag == 0);
  // shrink the size of w_inbag from n to n wts > 0
  this->w_inbag = w_inbag(rows_inbag);

 }

 void Tree::sample_cols(){

  // Start empty
  this->cols_node.set_size(mtry);
  uint cols_accepted = 0;

  // Set all to not selected
  std::vector<bool> temp;
  temp.resize(n_cols_total, false);

  std::uniform_int_distribution<uword> udist_cols(0, n_cols_total - 1);

  uword i, draw;

  for (i = 0; i < n_cols_total; ++i) {

   do {draw = udist_cols(random_number_generator); } while (temp[draw]);

   temp[draw] = true;

   if(is_col_splittable(draw)){
    cols_node[cols_accepted] = draw;
    cols_accepted++;
   }

   if(cols_accepted == mtry) break;

  }

  if(cols_accepted < mtry) cols_node.resize(cols_accepted);

 }

 void Tree::sample_cuts(){

  if(split_max_cuts >= cuts_all.size()){

   // no need for random sample if there are fewer valid cut-points
   // than the number of cut-points we planned to sample.
   cuts_sampled = cuts_all;

  } else { // split_max_cuts < cuts_all.size()

   cuts_sampled.resize(split_max_cuts);

   std::uniform_int_distribution<uword> udist_cuts(0, cuts_all.size() - 1);

   // Set all to not selected
   std::vector<bool> temp;
   temp.resize(cuts_all.size(), false);

   uword draw;

   for (uword i = 0; i < split_max_cuts; ++i) {

    do {draw = udist_cuts(random_number_generator); } while (temp[draw]);

    temp[draw] = true;

    cuts_sampled[i] = draw;

   }

   // important that cut-points are ordered from low to high
   cuts_sampled = cuts_all(cuts_sampled);
   cuts_sampled = sort(cuts_sampled);

  }

 }

 // not currently used but will be in the future
 // # nocov start
 bool Tree::is_col_splittable(uword j){

  uvec::iterator i;

  // initialize as 0 but do not make comparisons until x_first_value
  // is formally defined at the first instance of status == 1
  double x_first_value=0;

  bool x_first_undef = true;

  for (i = rows_node.begin(); i != rows_node.end(); ++i) {

   if(x_first_undef){

    x_first_value = x_inbag.at(*i, j);
    x_first_undef = false;

   } else {

    if(x_inbag.at(*i, j) != x_first_value){
     return(true);
    }

   }

  }

  if(verbosity > 4){
   // # nocov start
   mat x_print = x_inbag.rows(rows_node);
   Rcout << "   -- Column " << j << " was sampled but ";
   Rcout << "its unique values are " << unique(x_print.col(j));
   Rcout << std::endl;
   // # nocov end
  }

  return(false);

 }
 // # nocov end

 bool Tree::is_node_splittable(uword node_id){

  if(node_id == 0){

   // all inbag observations are in the first node
   rows_node = regspace<uvec>(0, n_rows_inbag-1);
   y_node = y_inbag;
   w_node = w_inbag;
   return(true);

  }

  rows_node = find(node_assignments == node_id);

  y_node = y_inbag.rows(rows_node);
  w_node = w_inbag(rows_node);

  bool result = is_node_splittable_internal();

  return(result);

 }

 // not currently used but will be in the future
 // # nocov start
 bool Tree::is_node_splittable_internal(){

  double n_obs = sum(w_node);

  return(n_obs >= 2*leaf_min_obs &&
         n_obs >= split_min_obs);

 }
 // # nocov end

 // not currently used but will be in the future
 // # nocov start
 void Tree::find_all_cuts(){

  // assume no valid cutpoints at first
  cuts_all.resize(0);

  uword i, j, k;

  uvec::iterator it, it_min, it_max;

  double n_obs = 0;

  // stop at end-1 b/c we access it+1 in lincomb_sort
  for(it = lincomb_sort.begin(); it < lincomb_sort.end()-1; ++it){

   n_obs += w_node[*it];

   // If we want to make the current value of lincomb a cut-point, we need
   // to make sure the next value of lincomb isn't equal to this current value.
   // Otherwise, we will have the same value of lincomb in both groups!

   if(lincomb[*it] != lincomb[*(it+1)]){

    if(n_obs >= leaf_min_obs) {

     if(verbosity > 3){
      // # nocov start
      Rcout << std::endl;
      Rcout << "   -- lower cutpoint: " << lincomb(*it) << std::endl;
      Rcout << "      - n_obs, left node:   " << n_obs   << std::endl;
      Rcout << std::endl;
      // # nocov end
     }

     break;

    }

   }

  }

  it_min = it;

  if(it == lincomb_sort.end()-1) {

   if(verbosity > 3){
    // # nocov start
    Rcout << "   -- Could not find a valid cut-point" << std::endl;
    // # nocov end
   }

   return;

  }

  // j = number of steps we have taken forward in lincomb
  j = it - lincomb_sort.begin();

  // reset before finding the upper limit
  n_obs=0;

  // stop at beginning+1 b/c we access it-1 in lincomb_sort
  for(it = lincomb_sort.end()-1; it >= lincomb_sort.begin()+1; --it){

   n_obs += w_node[*it];

   if(lincomb[*it] != lincomb[*(it-1)]){

    if(n_obs >= leaf_min_obs) {

     // the upper cutpoint needs to be one step below the current
     // it value, because we use x <= cp to determine whether a
     // value x goes to the left node versus the right node. So,
     // if it currently points to 3, and the next value down is 2,
     // then we want to say the cut-point is 2 because then all
     // values <= 2 will go left, and 3 will go right. This matters
     // when 3 is the highest value in the vector.

     --it;

     if(verbosity > 3){
      // # nocov start
      Rcout << std::endl;
      Rcout << "   -- upper cutpoint: " << lincomb(*it) << std::endl;
      Rcout << "      - n_obs, right node:   " << n_obs << std::endl;
      Rcout << std::endl;
      // # nocov end
     }

     break;

    }

   }

  }

  it_max = it;

  // k = n steps from beginning of sorted lincomb to current it
  k = it - lincomb_sort.begin();

  if(j > k){

   if(verbosity > 2) {
    // # nocov start
    Rcout << "   -- Could not find valid cut-points" << std::endl;
    // # nocov end
   }

   return;

  }

  // only one valid cutpoint
  if(j == k){ cuts_all = {j}; return; }

  i = 0;
  uvec output_middle(k-j);

  for(it = it_min+1;
      it < it_max; ++it){
   if(lincomb[*it] != lincomb[*(it+1)]){
    output_middle[i] = it - lincomb_sort.begin();
    i++;
   }
  }

  output_middle.resize(i);

  uvec output_left = {j};
  uvec output_right = {k};

  cuts_all = join_vert(output_left, output_middle, output_right);

 }
 // # nocov end

 double Tree::find_best_cut(){

  // initialize grouping for the current node
  // value of 1 indicates go to right node
  g_node.ones(lincomb.size());

  uvec::iterator it;

  uword it_start = 0, it_best = 0;

  double stat, stat_best = 0;

  if(verbosity > 3){
   // # nocov start
   Rcout << "   -- cutpoint (score)" << std::endl;
   // # nocov end
  }

  for(it = cuts_sampled.begin(); it != cuts_sampled.end(); ++it){

   // flip node assignments from left to right, up to the next cutpoint
   g_node.elem(lincomb_sort.subvec(it_start, *it)).fill(0);
   // compute split statistics with this cut-point
   stat = compute_split_score();
   // update leaderboard
   if(stat > stat_best) { stat_best = stat; it_best = *it; }
   // set up next loop run
   it_start = *it;

   if(verbosity > 3){
    // # nocov start
    Rcout << "   --- ";
    Rcout << lincomb.at(lincomb_sort(*it));
    Rcout << " (" << stat << "), ";
    Rcout << "N = " << sum(g_node % w_node) << " moving right";
    Rcout << std::endl;
    // # nocov end
   }

  }

  if(verbosity > 3){
   // # nocov start
   Rcout << std::endl;
   Rcout << "   -- best stat:  " << stat_best;
   Rcout << ", min to split: " << split_min_stat;
   Rcout << std::endl;
   Rcout << std::endl;
   // # nocov end
  }

  // do not split if best stat < minimum stat
  if(stat_best < split_min_stat){ return(R_PosInf); }

  // backtrack g_node to be what it was when best it was found
  if(it_best < it_start){
   g_node.elem(lincomb_sort.subvec(it_best+1, it_start)).fill(1);
  }

  // return the cut-point from best split
  return(lincomb[lincomb_sort[it_best]]);

 }

 void Tree::find_rows_inbag(arma::uword n_obs) {

  // it is assumed that:
  // - rows_oobag exists
  // - rows_oobag has length >= 1
  // - rows_oobag is sorted in ascending order

  rows_inbag.set_size(n_obs);
  uword rows_inbag_counter = 0;

  if(rows_oobag[0] != 0){
   for(arma::uword j = 0; j < rows_oobag.front(); ++j){
    rows_inbag[rows_inbag_counter] = j;
    rows_inbag_counter++;
   }
  }

  for(arma::uword i = 1; i < rows_oobag.size(); i++){
   if(rows_oobag[i-1]+1 != rows_oobag[i]){
    for(arma::uword j = rows_oobag[i-1]+1; j < rows_oobag[i]; ++j){
     rows_inbag[rows_inbag_counter] = j;
     rows_inbag_counter++;
    }
   }
  }

  if(rows_oobag.back() < n_obs){
   for(arma::uword j = rows_oobag.back()+1; j < n_obs; ++j){
    rows_inbag[rows_inbag_counter] = j;
    rows_inbag_counter++;
   }
  }

  rows_inbag.resize(rows_inbag_counter);

 }

 void Tree::sprout_leaf(uword node_id){

  if(verbosity > 2){
   // # nocov start
   Rcout << "-- sprouting node " << node_id << " into a leaf";
   Rcout << " (N = " << sum(w_node) << ")";
   Rcout << std::endl;
   Rcout << std::endl;
   // # nocov end
  }

  leaf_summary[node_id] = mean(y_node.col(0));

 }

 // not currently used but will be in the future
 // # nocov start
 double Tree::compute_split_score(){

  // default method is to pick one completely at random
  // (this won't stay the default - it's a placeholder)

  std::normal_distribution<double> ndist_score(0, 1);

  double result = ndist_score(random_number_generator);

  return(result);

 }
 // # nocov end

 // not currently used but will be in the future
 // # nocov start
 double Tree::compute_max_leaves(){

  // find maximum number of leaves for this tree
  // there are two ways to have maximal tree size:
  vec max_leaves_2ways = {
   //  1. every leaf node has exactly leaf_min_obs,
   n_obs_inbag / leaf_min_obs,
   //  2. every leaf node has exactly split_min_obs - 1,
   n_obs_inbag / (split_min_obs - 1)
  };

  double max_leaves = std::ceil(max(max_leaves_2ways));

  return(max_leaves);

 }
 // # nocov end

 void Tree::compute_oobag_vi(arma::vec* vi_numer,
                             VariableImportance vi_type) {

  resize_oobag_data();

  std::unique_ptr<Data> data_oobag { };
  data_oobag = std::make_unique<Data>(x_oobag, y_oobag, w_oobag);

  // using oobag = false for predict b/c data_oobag is already subsetted
  predict_leaf(data_oobag.get(), false);

  vec pred_values(data_oobag->n_rows);

  for(uword i = 0; i < pred_values.size(); ++i){
   pred_values[i] = leaf_summary[pred_leaf[i]];
  }

  // Compute normal prediction accuracy for each tree. Predictions already computed..
  double accuracy_normal = compute_prediction_accuracy(pred_values);

  if(verbosity > 1){
   // # nocov start
   Rcout << "  -- prediction accuracy before noising: ";
   Rcout << accuracy_normal << std::endl;
   Rcout << "  -- mean leaf pred: ";
   Rcout << mean(conv_to<vec>::from(pred_leaf));
   Rcout << std::endl << std::endl;
   // # nocov end
  }

  random_number_generator.seed(seed);

  // Randomly permute for all independent variables
  for (uword pred_col = 0; pred_col < data->get_n_cols(); ++pred_col) {

   // Check whether the i-th variable is used in the tree:
   bool pred_is_used = false;

   for(uint j = 0; j < coef_indices.size(); ++j){
    for(uword k = 0; k < coef_indices[j].size(); ++k){
     if(coef_indices[j][k] == pred_col){
      pred_is_used = true;
      break;
     }
    }
   }

   // proceed if the variable is used in the tree, otherwise vi = 0
   if (pred_is_used) {

    if(vi_type == VI_PERMUTE){
     // everyone gets the same permutation
     data_oobag->permute_col(pred_col, random_number_generator);
    } else if (vi_type == VI_NEGATE){
     negate_coef(pred_col);
    }

    predict_leaf(data_oobag.get(), false);

    for(uword i = 0; i < pred_values.size(); ++i){
     pred_values[i] = leaf_summary[pred_leaf[i]];
    }

    double accuracy_permuted = compute_prediction_accuracy(pred_values);

    if(verbosity > 3){
     // # nocov start
     Rcout << "   -- prediction accuracy after noising " << pred_col << ": ";
     Rcout << accuracy_permuted << std::endl;
     Rcout << "      - mean leaf pred: ";
     Rcout << mean(conv_to<vec>::from(pred_leaf));
     Rcout << std::endl << std::endl;
     // # nocov end
    }

    double accuracy_difference = accuracy_normal - accuracy_permuted;

    (*vi_numer)[pred_col] += accuracy_difference;

    if(vi_type == VI_PERMUTE){
     data_oobag->restore_col(pred_col);
    } else if (vi_type == VI_NEGATE){
     negate_coef(pred_col);
    }

   }
  }
 }

 void Tree::grow(arma::vec* vi_numer,
                 arma::uvec* vi_denom){

  this->vi_numer = vi_numer;
  this->vi_denom = vi_denom;

  sample_rows();

  // create inbag views of x, y, and w,
  this->x_inbag = data->x_rows(rows_inbag);
  this->y_inbag = data->y_rows(rows_inbag);

  this->n_obs_inbag = sum(w_inbag);
  this->n_rows_inbag = x_inbag.n_rows;

  node_assignments.zeros(n_rows_inbag);

  this->max_leaves = compute_max_leaves();
  this->max_nodes = (2 * max_leaves) - 1;

  if(verbosity > 2){
   // # nocov start
   Rcout << "- N obs inbag: " << n_obs_inbag;
   Rcout << std::endl;
   Rcout << "- N row inbag: " << n_rows_inbag;
   Rcout << std::endl;
   Rcout << "- max nodes: " << max_nodes;
   Rcout << std::endl;
   Rcout << "- max leaves: " << max_leaves;
   Rcout << std::endl;
   Rcout << std::endl;
   // # nocov end
  }

  // reserve memory for outputs (likely more than we need)
  cutpoint.resize(max_nodes);
  child_left.resize(max_nodes);
  coef_values.resize(max_nodes);
  coef_indices.resize(max_nodes);
  // memory for leaves based on corresponding tree type
  resize_leaves(max_nodes);

  // coordinate the order that nodes are grown.
  std::vector<uword> nodes_open;

  // start node 0
  nodes_open.push_back(0);

  // nodes to grow in the next run through the do-loop
  std::vector<uword> nodes_queued;

  // reserve space (most we could ever need is max_leaves)
  nodes_open.reserve(max_leaves);
  nodes_queued.reserve(max_leaves);

  // number of nodes in the tree
  uword n_nodes = 0;

  // iterate through nodes to be grown
  std::vector<uword>::iterator node;

  // ID of the left node (node_right = node_left + 1)
  uword node_left;

  do{

  for(node = nodes_open.begin(); node != nodes_open.end(); ++node){

   // determine rows in the current node and if it can be split
   if(!is_node_splittable(*node)){
    sprout_leaf(*node);
    continue;
   }

   uword n_retry = 0;

   // determines if a node is split or sprouted
   // (split means two new nodes are created)
   // (sprouted means the node becomes a leaf)
   for(; ;){

   // repeat until all the retries are spent.
    n_retry++;

    if(verbosity > 3){
     // # nocov start
     Rcout << "-- attempting to split node " << *node;
     Rcout << " (N = " << sum(w_node) << ",";
     Rcout << " try number " << n_retry << ")";
     Rcout << std::endl;
     Rcout << std::endl;
     // # nocov end
    }

    sample_cols();

    if(!cols_node.is_empty()){

     x_node = x_inbag(rows_node, cols_node);

     if(verbosity > 3) {
      // # nocov start
      print_uvec(cols_node, "columns sampled (showing up to 5)", 5);
      // # nocov end
     }

     // beta holds estimates (first item) and variance (second)
     // for the regression coefficients that created lincomb.
     // the variances are optional (only used for VI_ANOVA)
     mat beta;

     lincomb.zeros(x_node.n_rows);

     switch (lincomb_type) {

     case LC_NEWTON_RAPHSON: {

      beta = coxph_fit(x_node, y_node, w_node,
                       lincomb_scale, lincomb_ties_method,
                       lincomb_eps, lincomb_iter_max);

      break;

     }

     case LC_RANDOM_COEFS: {

      beta.set_size(x_node.n_cols, 1);

      std::uniform_real_distribution<double> unif_coef(0.0, 1.0);

      for(uword i = 0; i < x_node.n_cols; ++i){
       beta.at(i, 0) = unif_coef(random_number_generator);
      }

      break;

     }

     case LC_GLMNET: {

      NumericMatrix xx = wrap(x_node);
      NumericMatrix yy = wrap(y_node);
      NumericVector ww = wrap(w_node);

      // initialize function from tree object
      // (Functions can't be stored in C++ classes, but RObjects can)
      Function f_beta = as<Function>(lincomb_R_function);

      NumericMatrix beta_R = f_beta(xx, yy, ww,
                                    lincomb_alpha,
                                    lincomb_df_target);

      beta = mat(beta_R.begin(), beta_R.nrow(), beta_R.ncol(), false);

      break;

     }

     case LC_R_FUNCTION: {

      NumericMatrix xx = wrap(x_node);
      NumericMatrix yy = wrap(y_node);
      NumericVector ww = wrap(w_node);

      // initialize function from tree object
      // (Functions can't be stored in C++ classes, but RObjects can)
      Function f_beta = as<Function>(lincomb_R_function);

      NumericMatrix beta_R = f_beta(xx, yy, ww);

      beta = mat(beta_R.begin(), beta_R.nrow(), beta_R.ncol(), false);

      break;

     }

     } // end switch lincomb_type

     vec beta_est = beta.unsafe_col(0);

     if(verbosity > 3) {
      // # nocov start
      print_vec(beta_est, "linear combo weights (showing up to 5)", 5);
      // # nocov end
     }

     bool beta_all_zeros = find(beta_est != 0).is_empty();

     if(!beta_all_zeros){

      lincomb = x_node * beta_est;

      // sorted in ascending order
      lincomb_sort = sort_index(lincomb);

      // find all valid cutpoints for lincomb
      find_all_cuts();

      if(verbosity > 3 && cuts_all.is_empty()){
       // # nocov start
       Rcout << "   -- no cutpoints identified";
       Rcout << std::endl;
       // # nocov end
      }

      // empty cuts_all => no valid cutpoints => make leaf or retry
      if(!cuts_all.is_empty()){

       sample_cuts();

       double cut_point = find_best_cut();

       if(cut_point < R_PosInf){

        if(vi_type == VI_ANOVA && lincomb_type == LC_NEWTON_RAPHSON){

         // only do ANOVA variable importance when
         //  1. a split of the node is guaranteed
         //  2. the method used for lincombs allows it

         if(verbosity > 3){
          // # nocov start
          Rcout << "   -- p-values:" << std::endl;
          // # nocov end
         }

         vec beta_var = beta.unsafe_col(1);

         double pvalue;

         for(uword i = 0; i < beta_est.size(); ++i){

          (*vi_denom)[cols_node[i]]++;

          if(beta_est[i] != 0){

           pvalue = R::pchisq(pow(beta_est[i],2)/beta_var[i], 1, false, false);

           if(verbosity > 3){
            // # nocov start
            Rcout << "   --- column " << cols_node[i] << ": ";
            Rcout << pvalue;
            if(pvalue < 0.05) Rcout << "*";
            if(pvalue < 0.01) Rcout << "*";
            if(pvalue < 0.001) Rcout << "*";
            if(pvalue < vi_max_pvalue) Rcout << " [+1 to VI numerator]";
            Rcout << std::endl;
            // # nocov end
           }

           if(pvalue < vi_max_pvalue){ (*vi_numer)[cols_node[i]]++; }

          }

         }

         if(verbosity > 3){
          // # nocov start
          Rcout << std::endl;
          // # nocov end
         }

        }

        // make new nodes if a valid cutpoint was found
        node_left = n_nodes + 1;
        n_nodes += 2;
        // update tree parameters
        cutpoint[*node] = cut_point;
        coef_values[*node] = beta_est;
        coef_indices[*node] = cols_node;

        child_left[*node] = node_left;
        // re-assign observations in the current node
        // (note that g_node is 0 if left, 1 if right)
        node_assignments.elem(rows_node) = node_left + g_node;

        if(verbosity > 2){
         // # nocov start
         Rcout << "-- node " << *node << " was split into ";
         Rcout << "node " << node_left << " (left) and ";
         Rcout << node_left+1 << " (right)";
         Rcout << std::endl;
         Rcout << std::endl;
         // # nocov end
        }

        nodes_queued.push_back(node_left);
        nodes_queued.push_back(node_left + 1);
        break;

       }

      }

     }

    }

    if(n_retry >= split_max_retry){
     sprout_leaf(*node);
     break;
    }

   }


  }

  nodes_open = nodes_queued;
  nodes_queued.clear();

  } while (nodes_open.size() > 0);

  // don't forget to count the root node
  n_nodes++;

  cutpoint.resize(n_nodes);
  child_left.resize(n_nodes);
  coef_values.resize(n_nodes);
  coef_indices.resize(n_nodes);

  resize_leaves(n_nodes);

 } // Tree::grow

 void Tree::predict_leaf(Data* prediction_data, bool oobag) {

  pred_leaf.zeros(prediction_data->n_rows);

  // if tree is root node, 0 is the correct leaf prediction
  if(coef_values.size() == 0) return;

  if(verbosity > 2){
   // # nocov start
   Rcout << "   -- computing leaf predictions" << std::endl;
   // # nocov end
  }

  uvec obs_in_node;

  // it iterates over the observations in a node
  uvec::iterator it;

  // i iterates over nodes, j over observations
  // uword i, j;

  for(uword i = 0; i < coef_values.size(); i++){

   // if child_left == 0, it's a leaf (no need to find next child)
   if(child_left[i] != 0){

    if(i == 0 && oobag){
     obs_in_node = rows_oobag;
    } else if (i == 0 && !oobag) {
     obs_in_node = regspace<uvec>(0, 1, pred_leaf.size()-1);
    } else {
     obs_in_node = find(pred_leaf == i);
    }

    if(obs_in_node.size() > 0){

     x_node = prediction_data->x_submat(obs_in_node, coef_indices[i]);

     it = obs_in_node.begin();

     for(uword j = 0; j < obs_in_node.size(); ++j, ++it){

      if(dot(x_node.row(j), coef_values[i]) <= cutpoint[i]) {

       pred_leaf[*it] = child_left[i];

      } else {

       pred_leaf[*it] = child_left[i]+1;

      }

     }

     if(verbosity > 4){
      // # nocov start
      uvec in_left = find(pred_leaf == child_left[i]);
      uvec in_right = find(pred_leaf == child_left[i]+1);
      Rcout << "No. to node " << child_left[i] << ": ";
      Rcout << in_left.size() << "; " << std::endl;
      Rcout << "No. to node " << child_left[i]+1 << ": ";
      Rcout << in_right.size() << std::endl << std::endl;
      // # nocov end
     }

    }

   }

  }

  if(oobag){ pred_leaf.elem(rows_inbag).fill(max_nodes); }

 }

 double Tree::compute_prediction_accuracy(arma::vec& preds){

  if (oobag_eval_type == EVAL_R_FUNCTION){

   NumericMatrix y_wrap = wrap(y_oobag);
   NumericVector w_wrap = wrap(w_oobag);
   NumericVector p_wrap = wrap(preds);

   // initialize function from tree object
   // (Functions can't be stored in C++ classes, but RObjects can)
   Function f_oobag = as<Function>(oobag_R_function);

   NumericVector result_R = f_oobag(y_wrap, w_wrap, p_wrap);

   return(result_R[0]);

  }

  return(compute_prediction_accuracy_internal(preds));

 }

 void Tree::negate_coef(arma::uword pred_col){

  for(uint j = 0; j < coef_indices.size(); ++j){

   for(uword k = 0; k < coef_indices[j].size(); ++k){
    if(coef_indices[j][k] == pred_col){
     coef_values[j][k] *= (-1);
    }
   }

  }

 }

 } // namespace aorsf

