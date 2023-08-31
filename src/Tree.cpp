/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "Tree.h"
#include "Coxph.h"
#include "NodeSplitStats.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 Tree::Tree() :
   data(0),
   n_cols_total(0),
   n_rows_total(0),
   seed(0),
   mtry(0),
   vi_type(VI_NONE),
   vi_max_pvalue(DEFAULT_ANOVA_VI_PVALUE),
   leaf_min_events(DEFAULT_LEAF_MIN_EVENTS),
   leaf_min_obs(DEFAULT_LEAF_MIN_OBS),
   split_rule(DEFAULT_SPLITRULE),
   split_min_events(DEFAULT_SPLIT_MIN_EVENTS),
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
   lincomb_R_function(0) {

 }

 Tree::Tree(std::vector<double>& cutpoint,
            std::vector<arma::uword>& child_left,
            std::vector<arma::vec>& coef_values,
            std::vector<arma::uvec>& coef_indices,
            std::vector<arma::vec>& leaf_pred_horizon,
            std::vector<arma::vec>& leaf_pred_surv,
            std::vector<arma::vec>& leaf_pred_chf) :
 data(0),
 n_cols_total(0),
 n_rows_total(0),
 seed(0),
 mtry(0),
 vi_type(VI_NONE),
 vi_max_pvalue(DEFAULT_ANOVA_VI_PVALUE),
 leaf_min_events(DEFAULT_LEAF_MIN_EVENTS),
 leaf_min_obs(DEFAULT_LEAF_MIN_OBS),
 split_rule(DEFAULT_SPLITRULE),
 split_min_events(DEFAULT_SPLIT_MIN_EVENTS),
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
 cutpoint(cutpoint),
 child_left(child_left),
 coef_values(coef_values),
 coef_indices(coef_indices),
 leaf_pred_horizon(leaf_pred_horizon),
 leaf_pred_surv(leaf_pred_surv),
 leaf_pred_chf(leaf_pred_chf) {

 }


 void Tree::init(Data* data,
                 int seed,
                 arma::uword mtry,
                 double leaf_min_events,
                 double leaf_min_obs,
                 VariableImportance vi_type,
                 double vi_max_pvalue,
                 SplitRule split_rule,
                 double split_min_events,
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
                 RObject lincomb_R_function){

  // Initialize random number generator and set seed
  random_number_generator.seed(seed);

  this->data = data;
  this->n_cols_total = data->n_cols;
  this->n_rows_total = data->n_rows;
  this->seed = seed;
  this->mtry = mtry;
  this->leaf_min_events = leaf_min_events;
  this->leaf_min_obs = leaf_min_obs;
  this->vi_type = vi_type;
  this->vi_max_pvalue = vi_max_pvalue;
  this->split_rule = split_rule;
  this->split_min_events = split_min_events;
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

 }

 void Tree::sample_rows(){

  uword i, draw, n = data->n_rows;

  // Start with all samples OOB
  vec w_inbag(n, fill::zeros);

  std::uniform_int_distribution<uword> unif_dist(0, n - 1);

  // sample with replacement
  for (i = 0; i < n; ++i) {
   draw = unif_dist(random_number_generator);
   ++w_inbag[draw];
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
  std::vector<uword> cols_accepted;
  cols_accepted.reserve(mtry);

  // Set all to not selected
  std::vector<bool> temp;
  temp.resize(n_cols_total, false);

  std::uniform_int_distribution<uword> unif_dist(0, n_cols_total - 1);

  uword i, draw;

  for (i = 0; i < n_cols_total; ++i) {

   do { draw = unif_dist(random_number_generator); } while (temp[draw]);

   temp[draw] = true;

   if(is_col_splittable(draw)){
    cols_accepted.push_back(draw);
   }

   if(cols_accepted.size() == mtry) break;

  }

  this->cols_node = uvec(cols_accepted.data(),
                         cols_accepted.size(),
                         false,
                         true);

 }

 bool Tree::is_col_splittable(uword j){

  uvec::iterator i;

  // initialize as 0 but do not make comparisons until x_first_value
  // is formally defined at the first instance of status == 1
  double x_first_value=0;

  bool x_first_undef = true;

  for (i = rows_node.begin(); i != rows_node.end(); ++i) {

   // if event occurred for this observation
   if(y_inbag.at(*i, 1) == 1){

    if(x_first_undef){

     x_first_value = x_inbag.at(*i, j);
     x_first_undef = false;

    } else {

     if(x_inbag.at(*i, j) != x_first_value){
      return(true);
     }

    }

   }

  }

  if(VERBOSITY > 1){

   mat x_print = x_inbag.rows(rows_node);
   mat y_print = y_inbag.rows(rows_node);

   uvec rows_event = find(y_print.col(1) == 1);
   x_print = x_print.rows(rows_event);

   Rcout << "Column " << j << " was sampled but ";
   Rcout << "unique values of column " << j << " are ";
   Rcout << unique(x_print.col(j)) << std::endl;

  }

  return(false);

 }

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

  double n_risk = sum(w_node);
  double n_events = sum(y_node.col(1) % w_node);

  return(n_events >= 2*leaf_min_events &&
         n_risk   >= 2*leaf_min_obs &&
         n_events >=   split_min_events &&
         n_risk   >=   split_min_obs);

 }


 uvec Tree::find_cutpoints(){

  vec y_status = y_node.unsafe_col(1);

  // placeholder with values indicating invalid cps
  uvec output;

  uword i, j, k;

  uvec::iterator it, it_min, it_max;

  double n_events = 0, n_risk = 0;

  if(VERBOSITY > 1){
   Rcout << "----- finding lower bound for cut-points -----" << std::endl;
  }

  // stop at end-1 b/c we access it+1 in lincomb_sort
  for(it = lincomb_sort.begin(); it < lincomb_sort.end()-1; ++it){

   n_events += y_status[*it] * w_node[*it];
   n_risk += w_node[*it];


   if(VERBOSITY > 2){
    Rcout << "current value: "<< lincomb(*it) << " -- ";
    Rcout << "next: "<< lincomb(*(it+1))      << " -- ";
    Rcout << "N events: " << n_events         << " -- ";
    Rcout << "N risk: "   << n_risk           << std::endl;
   }

   // If we want to make the current value of lincomb a cut-point, we need
   // to make sure the next value of lincomb isn't equal to this current value.
   // Otherwise, we will have the same value of lincomb in both groups!

   if(lincomb[*it] != lincomb[*(it+1)]){

    if( n_events >= leaf_min_events &&
        n_risk   >= leaf_min_obs) {

     if(VERBOSITY > 0){
      Rcout << std::endl;
      Rcout << "lower cutpoint: "         << lincomb(*it) << std::endl;
      Rcout << " - n_events, left node: " << n_events << std::endl;
      Rcout << " - n_risk, left node:   " << n_risk   << std::endl;
      Rcout << std::endl;
     }

     break;

    }

   }

  }

  it_min = it;

  if(it == lincomb_sort.end()-1) {

   if(VERBOSITY > 1){
    Rcout << "Could not find a valid cut-point" << std::endl;
   }

   return(output);

  }

  // j = number of steps we have taken forward in lincomb
  j = it - lincomb_sort.begin();

  // reset before finding the upper limit
  n_events=0, n_risk=0;

  if(VERBOSITY > 1){
   Rcout << "----- finding upper bound for cut-points -----" << std::endl;
  }

  // stop at beginning+1 b/c we access it-1 in lincomb_sort
  for(it = lincomb_sort.end()-1; it >= lincomb_sort.begin()+1; --it){

   n_events += y_status[*it] * w_node[*it];
   n_risk   += w_node[*it];

   if(VERBOSITY > 2){
    Rcout << "current value: "<< lincomb(*it)  << " ---- ";
    Rcout << "next value: "<< lincomb(*(it-1)) << " ---- ";
    Rcout << "N events: " << n_events       << " ---- ";
    Rcout << "N risk: " << n_risk           << std::endl;
   }

   if(lincomb(*it) != lincomb(*(it-1))){

    if( n_events >= leaf_min_events &&
        n_risk   >= leaf_min_obs ) {

     // the upper cutpoint needs to be one step below the current
     // it value, because we use x <= cp to determine whether a
     // value x goes to the left node versus the right node. So,
     // if it currently points to 3, and the next value down is 2,
     // then we want to say the cut-point is 2 because then all
     // values <= 2 will go left, and 3 will go right. This matters
     // when 3 is the highest value in the vector.

     --it;

     if(VERBOSITY > 0){
      Rcout << std::endl;
      Rcout << "upper cutpoint: " << lincomb(*it) << std::endl;
      Rcout << " - n_events, right node: " << n_events    << std::endl;
      Rcout << " - n_risk, right node:   " << n_risk      << std::endl;
      Rcout << std::endl;
     }

     break;

    }

   }

  }

  it_max = it;

  // k = n steps from beginning of sorted lincomb to current it
  k = it - lincomb_sort.begin();

  if(j > k){

   if(VERBOSITY > 0) {
    Rcout << "Could not find valid cut-points" << std::endl;
   }

   return(output);

  }

  // only one valid cutpoint
  if(j == k){

   uvec output = {j};
   return(output);

  }

  i = 0;
  uvec output_middle(k-j);

  for(it = it_min+1;
      it < it_max; ++it){
   if(lincomb(*it) != lincomb(*(it+1))){
    output_middle[i] = it - lincomb_sort.begin();
    i++;
   }
  }

  output_middle.resize(i);

  uvec output_left = {j};
  uvec output_right = {k};

  output = join_vert(output_left, output_middle, output_right);

  return(output);

 }

 double Tree::score_logrank(){

  double
  n_risk=0,
   g_risk=0,
   observed=0,
   expected=0,
   V=0,
   temp1,
   temp2,
   n_events;

  vec y_time = y_node.unsafe_col(0);
  vec y_status = y_node.unsafe_col(1);

  bool break_loop = false;

  uword i = y_node.n_rows-1;

  // breaking condition of outer loop governed by inner loop
  for (; ;){

   temp1 = y_time[i];

   n_events = 0;

   for ( ; y_time[i] == temp1; i--) {

    n_risk += w_node[i];
    n_events += y_status[i] * w_node[i];
    g_risk += g_node[i] * w_node[i];
    observed += y_status[i] * g_node[i] * w_node[i];

    if(i == 0){
     break_loop = true;
     break;
    }

   }

   // should only do these calculations if n_events > 0,
   // but in practice its often faster to multiply by 0
   // versus check if n_events is > 0.

   temp2 = g_risk / n_risk;
   expected += n_events * temp2;

   // update variance if n_risk > 1 (if n_risk == 1, variance is 0)
   // definitely check if n_risk is > 1 b/c otherwise divide by 0
   if (n_risk > 1){
    temp1 = n_events * temp2 * (n_risk-n_events) / (n_risk-1);
    V += temp1 * (1 - temp2);
   }

   if(break_loop) break;

  }

  return(pow(expected-observed, 2) / V);

 }

 double Tree::node_split(arma::uvec& cuts_all){

  // sample a subset of cutpoints.
  uvec cuts_sampled;

  if(split_max_cuts >= cuts_all.size()){

   // no need for random sample if there are fewer valid cut-points
   // than the number of cut-points we planned to sample.
   cuts_sampled = cuts_all;

  } else { // split_max_cuts < cuts_all.size()

   // cuts_sampled.set_size(split_max_cuts);
   //
   // std::uniform_int_distribution<uword> unif_dist(0, cuts_all.size() - 1);
   //
   // // sample without replacement
   // for (uword i = 0; i < split_max_cuts; ++i) {
   //
   //  uword draw = unif_dist(random_number_generator);
   //
   //  // Ensure the drawn number is not already in the sample
   //  while (std::find(cuts_sampled.begin(),
   //                   cuts_sampled.end(),
   //                   cuts_all[draw]) != cuts_sampled.end()) {
   //
   //   draw = unif_dist(random_number_generator);
   //
   //  }
   //
   //  cuts_sampled[i] = cuts_all[draw];
   //
   // }
   //
   //
   //
   // // important that cut-points are ordered from low to high
   // cuts_sampled = sort(cuts_sampled);
   //
   // if(VERBOSITY > 1){
   //
   //  Rcout << "Randomly sampled cutpoints: ";
   //  Rcout << std::endl;
   //  Rcout << lincomb(lincomb_sort(cuts_sampled));
   //  Rcout << std::endl;
   //  Rcout << std::endl;
   //
   // }

   // non-random version
   cuts_sampled = linspace<uvec>(cuts_all.front(),
                                 cuts_all.back(),
                                 split_max_cuts);

  }


  // initialize grouping for the current node
  // value of 1 indicates go to right node
  g_node.ones(lincomb.size());

  uvec::iterator it;

  uword it_start = 0, it_best;

  double stat, stat_best = 0;

  for(it = cuts_sampled.begin(); it != cuts_sampled.end(); ++it){

   // flip node assignments from left to right, up to the next cutpoint
   g_node.elem(lincomb_sort.subvec(it_start, *it)).fill(0);
   // compute split statistics with this cut-point
   stat = score_logrank();
   // update leaderboard
   if(stat > stat_best) { stat_best = stat; it_best = *it; }
   // set up next loop run
   it_start = *it;

   if(VERBOSITY > 1){
    mat temp = join_rows(lincomb, conv_to<vec>::from(g_node));
    temp = join_rows(temp, w_node);
    temp = temp.rows(lincomb_sort);
    Rcout << "testing cutpoint: " << lincomb.at(lincomb_sort(*it));
    Rcout << std::endl;
    Rcout << std::endl;
    print_mat(temp, "lincomb & g_node & w_node", 20, 20);
    Rcout << "logrank stat for this cutpoint: " << stat;
    Rcout << std::endl;
    Rcout << " ------------------------------------------- ";
    Rcout << std::endl;
    Rcout << std::endl;
   }

  }

  // do not split if best stat < minimum stat
  if(stat_best < split_min_stat){

   if(VERBOSITY > 1){
    Rcout << "best split stat, " << stat_best;
    Rcout << ", was < split_min_stat, " << split_min_stat;
    Rcout << std::endl;
   }

   return(R_PosInf);

  }

  // backtrack g_node to be what it was when best it was found
  if(it_best < it_start){
   g_node.elem(lincomb_sort.subvec(it_best+1, it_start)).fill(1);
  }


  // return the cut-point from best split
  return(lincomb[lincomb_sort[it_best]]);

 }

 void Tree::node_sprout(uword node_id){

  if(VERBOSITY > 0){
   Rcout << "sprouting new leaf with node " << node_id;
   Rcout << std::endl;
   Rcout << std::endl;
  }

  // reserve as much size as could be needed (probably more)
  mat leaf_data(y_node.n_rows, 3);

  uword person = 0;

  // find the first unique event time
  while(y_node.at(person, 1) == 0 && person < y_node.n_rows){
   person++;
  }

  // person corresponds to first event or last censor time
  leaf_data.at(0, 0) = y_node.at(person, 0);

  // if no events in this node:
  if(person == y_node.n_rows){

   vec temp_surv(1, arma::fill::ones);
   vec temp_chf(1, arma::fill::zeros);

   leaf_pred_horizon[node_id] = leaf_data.col(0);
   leaf_pred_surv[node_id] = temp_surv;
   leaf_pred_chf[node_id] = temp_chf;

   return;

  }

  double temp_time = y_node.at(person, 0);

  uword i = 1;

  // find the rest of the unique event times
  for( ; person < y_node.n_rows; person++){

   if(temp_time != y_node.at(person, 0) && y_node.at(person, 1) == 1){

    leaf_data.at(i, 0) = y_node.at(person,0);
    temp_time = y_node.at(person, 0);
    i++;

   }

  }

  leaf_data.resize(i, 3);

  // reset for kaplan meier loop
  person = 0; i = 0;
  double n_risk = sum(w_node);
  double temp_surv = 1.0;
  double temp_haz = 0.0;

  do {

   double n_events   = 0;
   double n_risk_sub = 0;
   temp_time = y_node.at(person, 0);

   while(y_node.at(person, 0) == temp_time){

    n_risk_sub += w_node.at(person);
    n_events += y_node.at(person, 1) * w_node.at(person);

    if(person == y_node.n_rows-1) break;

    person++;

   }

   // only do km if a death was observed

   if(n_events > 0){

    temp_surv = temp_surv * (n_risk - n_events) / n_risk;

    temp_haz = temp_haz + n_events / n_risk;

    leaf_data.at(i, 1) = temp_surv;
    leaf_data.at(i, 2) = temp_haz;
    i++;

   }

   n_risk -= n_risk_sub;

  } while (i < leaf_data.n_rows);


  if(VERBOSITY > 1) print_mat(leaf_data, "leaf_data", 10, 5);

  leaf_pred_horizon[node_id] = leaf_data.col(0);
  leaf_pred_surv[node_id] = leaf_data.col(1);
  leaf_pred_chf[node_id] = leaf_data.col(2);

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
  this->n_events_inbag = sum(w_inbag % y_inbag.col(1));
  this->n_rows_inbag = x_inbag.n_rows;

  if(VERBOSITY > 0){

   Rcout << "Effective sample size: " << n_obs_inbag;
   Rcout << std::endl;
   Rcout << "Effective number of events: " << n_events_inbag;
   Rcout << std::endl;
   Rcout << "Number of unique rows in x: " << n_rows_inbag;
   Rcout << std::endl;
   Rcout << std::endl;

  }

  node_assignments.zeros(n_rows_inbag);

  // find maximum number of leaves for this tree
  // there are four ways to have maximal tree size:
  vec max_leaves_4ways = {
   //  1. every leaf node has exactly leaf_min_obs,
   n_obs_inbag / leaf_min_obs,
   //  2. every leaf node has exactly leaf_min_events,
   n_events_inbag / leaf_min_events,
   //  3. every leaf node has exactly split_min_obs - 1,
   n_obs_inbag / (split_min_obs - 1),
   //  4. every leaf node has exactly split_min_events-1
   n_events_inbag / (split_min_events - 1)
  };

  // number of nodes total in binary tree is 2*L - 1,
  // where L is the number of leaf nodes in the tree.
  // (can prove by induction)
  double max_leaves = std::ceil(max(max_leaves_4ways));
  double max_nodes = (2 * max_leaves) - 1;

  if(VERBOSITY > 0){

   Rcout << "Max number of nodes for this tree: " << max_nodes;
   Rcout << std::endl;
   Rcout << "Max number of leaves for this tree: " << max_leaves;
   Rcout << std::endl;
   Rcout << std::endl;


  }

  // reserve memory for outputs (likely more than we need)
  cutpoint.resize(max_nodes);
  child_left.resize(max_nodes);
  coef_values.resize(max_nodes);
  coef_indices.resize(max_nodes);
  leaf_pred_horizon.resize(max_nodes);
  leaf_pred_surv.resize(max_nodes);
  leaf_pred_chf.resize(max_nodes);

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

  // all possible cut-points for a linear combination
  uvec cuts_all;

  do{

  for(node = nodes_open.begin(); node != nodes_open.end(); ++node){

   if(VERBOSITY > 0){
    Rcout << "growing node " << *node;
    Rcout << std::endl << std::endl;
   }


   // determine rows in the current node and if it can be split
   if(!is_node_splittable(*node)){

    node_sprout(*node);
    continue;

   }

   uword n_retry = 0;

   // determines if a node is split or sprouted
   // (split means two new nodes are created)
   // (sprouted means the node becomes a leaf)
   for(; ;){

   // repeat until all the retries are spent.
    n_retry++;

    if(VERBOSITY > 1){

     Rcout << "beginning try no. " << n_retry;
     Rcout << std::endl << std::endl;

    }

    sample_cols();

    x_node = x_inbag(rows_node, cols_node);

    if(VERBOSITY > 1) {
     print_mat(x_node, "x_node", 20, 20);
     print_mat(y_node, "y_node", 20, 20);
    }

    // beta holds estimates (first item) and variance (second)
    // for the regression coefficients that created lincomb.
    // the variances are optional (only used for VI_ANOVA)
    mat beta;

    lincomb.zeros(x_node.n_rows);

    switch (lincomb_type) {

    case NEWTON_RAPHSON: {

     beta = coxph_fit(x_node, y_node, w_node,
                      lincomb_scale, lincomb_ties_method,
                      lincomb_eps, lincomb_iter_max);


     break;

    }

    case RANDOM_COEFS: {

     beta.set_size(x_node.n_cols, 1);

     std::uniform_real_distribution<double> unif_coef(0.0, 1.0);

     for(uword i = 0; i < x_node.n_cols; ++i){
      beta.at(i, 0) = unif_coef(random_number_generator);
     }

     break;

    }

    case R_FUNCTION: {

     // NumericMatrix xx = ;
     // NumericMatrix yy = ;
     // NumericVector ww = ;

     // initialize function from tree object
     // (Functions can't be stored in C++ classes, but RObjects can)
     Function f_beta = as<Function>(lincomb_R_function);

     NumericMatrix beta_R = f_beta(wrap(x_node),
                                   wrap(y_node),
                                   wrap(w_node));

     beta = mat(beta_R.begin(), beta_R.nrow(), beta_R.ncol(), false);

     break;

    }

    } // switch lincomb_type

    vec beta_est = beta.unsafe_col(0);

    lincomb = x_node * beta_est;

    // sorted in ascending order
    lincomb_sort = sort_index(lincomb);

    // find all valid cutpoints for lincomb
    cuts_all = find_cutpoints();

    // empty cuts_all => no valid cutpoints => make leaf or retry
    if(!cuts_all.is_empty()){

     double cut_point = node_split(cuts_all);

     if(cut_point < R_PosInf){

      if(vi_type == VI_ANOVA && lincomb_type == NEWTON_RAPHSON){

       // only do ANOVA variable importance when
       //  1. a split of the node is guaranteed
       //  2. the method used for lincombs allows it

       vec beta_var = beta.unsafe_col(1);

       double pvalue;

       for(uword i = 0; i < beta_est.size(); ++i){

        (*vi_denom)[cols_node[i]]++;

        if(beta_est[i] != 0){

         pvalue = R::pchisq(pow(beta_est[i],2)/beta_var[i], 1, false, false);

         if(pvalue < vi_max_pvalue){ (*vi_numer)[cols_node[i]]++; }

        }

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

      if(VERBOSITY > 1){
       Rcout << "Split successful: unique node assignments: ";
       Rcout << std::endl;
       Rcout << unique(node_assignments).t();
       Rcout << std::endl;
      }

      nodes_queued.push_back(node_left);
      nodes_queued.push_back(node_left + 1);
      break;

     }

    }

    if(n_retry == split_max_retry){
     node_sprout(*node);
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
  leaf_pred_horizon.resize(n_nodes);
  leaf_pred_surv.resize(n_nodes);
  leaf_pred_chf.resize(n_nodes);

 } // Tree::grow

 void Tree::predict_leaf(Data* prediction_data, bool oobag) {

  pred_leaf.zeros(prediction_data->n_rows);

  if(VERBOSITY > 0){
   Rcout << "---- computing leaf predictions ----" << std::endl;
  }

  uvec obs_in_node;
  // it iterates over the observations in a node
  uvec::iterator it;

  // i iterates over nodes, j over observations
  uword i, j;

  for(i = 0; i < coef_values.size(); i++){

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

     lincomb = prediction_data->x_submat(obs_in_node, coef_indices[i]) * coef_values[i];

     it = obs_in_node.begin();

     for(j = 0; j < lincomb.size(); ++j, ++it){

      if(lincomb[j] <= cutpoint[i]) {

       pred_leaf[*it] = child_left[i];

      } else {

       pred_leaf[*it] = child_left[i]+1;

      }

     }

     if(VERBOSITY > 0){

      uvec in_left = find(pred_leaf == child_left[i]);
      uvec in_right = find(pred_leaf == child_left[i]+1);

      Rcout << "No. to node " << child_left[i] << ": ";
      Rcout << in_left.size() << "; " << std::endl;
      Rcout << "No. to node " << child_left[i]+1 << ": ";
      Rcout << in_right.size() << std::endl << std::endl;

     }

    }

   }

  }

 }

 void Tree::predict_value(arma::mat* pred_output,
                          arma::vec* pred_denom,
                          arma::vec& pred_times,
                          char pred_type,
                          bool oobag){

  uvec pred_leaf_sort = sort_index(pred_leaf, "ascend");

  uvec::iterator it = pred_leaf_sort.begin();

  // oobag leaf prediction has zeros for inbag rows
  if(oobag){ while(pred_leaf[*it] == 0){ ++it; } }

  double pred_t0;

  if(pred_type == 'S' || pred_type == 'R'){
   pred_t0 = 1;
  } else {
   pred_t0 = 0;
  }

  uword i, j;

  vec leaf_times, leaf_values;

  vec temp_vec(pred_times.size());
  double temp_dbl;

  do {

   uword leaf_id = pred_leaf[*it];

   // copies of leaf data using same aux memory
   leaf_times = vec(leaf_pred_horizon[leaf_id].begin(),
                    leaf_pred_horizon[leaf_id].size(),
                    false);

   leaf_values = vec(leaf_pred_surv[leaf_id].begin(),
                     leaf_pred_surv[leaf_id].size(),
                     false);

   if(leaf_values.is_empty()) Rcpp::stop("empty leaf");

   // don't reset i in the loop.
   // (wasteful b/c leaf_times ascend)
   i = 0;

   for(j = 0; j < pred_times.size(); j++){

    // t is the current prediction time
    double t = pred_times[j];

   // if t < t', where t' is the max time in this leaf,
   // then we may find a time t* such that t* < t < t'.
   // If so, prediction should be anchored to t*.
   // But, there may be multiple t* < t, and we want to
   // find the largest t* that is < t, so we find the
   // first t** > t and assign t* to be whatever came
   // right before t**.
   if(t < leaf_times.back()){

    for(; i < leaf_times.size(); i++){

     // we found t**
     if (leaf_times[i] > t){

      if(i == 0)
       // first leaf event occurred after prediction time
       temp_dbl = pred_t0;
      else
       // t* is the time value just before t**, so use i-1
       temp_dbl = leaf_values[i-1];

      break;

     } else if (leaf_times[i] == t){
      // pred_horizon just happens to equal a leaf time
      temp_dbl = leaf_values[i];

      break;

     }

    }

   } else {
    // if t > t' use the last recorded prediction
    temp_dbl = leaf_values.back();

   }

   temp_vec[j] = temp_dbl;

   }

   (*pred_output).row(*it) += temp_vec.t();
   if(oobag) (*pred_denom)[*it]++;
   ++it;

   if(it < pred_leaf_sort.end()){

    while(leaf_id == pred_leaf[*it]){

     (*pred_output).row(*it) += temp_vec.t();
     if(oobag) (*pred_denom)[*it]++;

     if (it == pred_leaf_sort.end()-1){ break; } else { ++it; }

    }

   }

  } while (it < pred_leaf_sort.end());



 }



 } // namespace aorsf

