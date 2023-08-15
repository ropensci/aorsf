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

 Tree::Tree(){ }

 void Tree::init(Data* data,
                 int seed,
                 arma::uword mtry,
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
                 bool lincomb_scale,
                 double lincomb_alpha,
                 arma::uword lincomb_df_target,
                 arma::uword lincomb_ties_method){


  this->data = data;
  this->n_cols_total = data->n_cols;

  this->seed = seed;
  this->mtry = mtry;
  this->leaf_min_events = leaf_min_events;
  this->leaf_min_obs = leaf_min_obs;
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

 }

 void Tree::sample_rows(){

  // Initialize random number generator and set seed
  random_number_generator.seed(seed);

  uword i, draw, n = data->n_rows;

  // Start with all samples OOB
  vec boot_wts(n, fill::zeros);

  std::uniform_int_distribution<uword> unif_dist(0, n - 1);

  // sample with replacement
  for (i = 0; i < n; ++i) {
   draw = unif_dist(random_number_generator);
   ++boot_wts[draw];
  }

  // multiply boot_wts by user specified weights.
  if(data->has_weights){
   boot_wts = boot_wts = boot_wts % data->w;
  }

  uvec rows_inbag = find(boot_wts > 0);
  this->x_inbag = data->x_rows(rows_inbag);
  this->y_inbag = data->y_rows(rows_inbag);
  this->w_inbag = data->w_subvec(rows_inbag);
  this->rows_oobag = find(boot_wts == 0);

  if(VERBOSITY > 0){

   print_mat(x_inbag, "x_inbag", 5, 5);

  }

  rows_inbag.clear();
  boot_wts.clear();

 }

 void Tree::sample_cols(){

  // Start empty
  std::vector<uword> cols_assessed, cols_accepted;
  cols_assessed.reserve(n_cols_total);
  cols_accepted.reserve(mtry);

  std::uniform_int_distribution<uword> unif_dist(0, n_cols_total - 1);

  uword i, draw;

  for (i = 0; i < n_cols_total; ++i) {

   draw = unif_dist(random_number_generator);

   // Ensure the drawn number is not already in the sample
   while (std::find(cols_assessed.begin(),
                    cols_assessed.end(),
                    draw) != cols_assessed.end()) {

    draw = unif_dist(random_number_generator);


   }

   cols_assessed.push_back(draw);

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

  vec status = y_inbag.unsafe_col(1);

  uvec::iterator i;

  // initialize as 0 but do not make comparisons until x_first_value
  // is formally defined at the first instance of status == 1
  double x_first_value=0;

  bool x_first_undef = true;

  for (i = rows_node.begin(); i != rows_node.end(); ++i) {

   if(status[*i] == 1){

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

  return(false);

 }

 bool Tree::is_left_node_splittable(){

  return(n_events_left >= 2*leaf_min_events &&
         n_risk_left   >= 2*leaf_min_obs &&
         n_events_left >=   split_min_events &&
         n_risk_left   >=   split_min_obs);

 }

 bool Tree::is_right_node_splittable(){

  return(n_events_right >= 2*leaf_min_events &&
         n_risk_right   >= 2*leaf_min_obs &&
         n_events_right >=   split_min_events &&
         n_risk_right   >=   split_min_obs);

 }

 uvec Tree::find_cutpoints(vec& lincomb, uvec& lincomb_sort){

  vec y_status = y_node.unsafe_col(1);

  // placeholder with values indicating invalid cps
  uvec output;

  uword i, j, k;

  uvec::iterator it, it_min, it_max;

  double n_events = 0, n_risk = 0;

  if(VERBOSITY > 0){
   Rcout << "----- finding lower bound for cut-points -----" << std::endl;
  }

  // stop at end-1 b/c we access it+1 in lincomb_sort
  for(it = lincomb_sort.begin(); it < lincomb_sort.end()-1; ++it){

   n_events += y_status[*it] * w_node[*it];
   n_risk += w_node[*it];


   if(VERBOSITY > 2){
    Rcout << "current value: "<< lincomb(*it)  << " ---- ";
    Rcout << "next value: "<< lincomb(*(it+1)) << " ---- ";
    Rcout << "N events: " << n_events       << " ---- ";
    Rcout << "N risk: " << n_risk           << std::endl;
   }

   // If we want to make the current value of lincomb a cut-point, we need
   // to make sure the next value of lincomb isn't equal to this current value.
   // Otherwise, we will have the same value of lincomb in both groups!

   if(lincomb[*it] != lincomb[*(it+1)]){

    if( n_events >= leaf_min_events &&
        n_risk   >= leaf_min_obs) {

     if(VERBOSITY > 1){
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

  if(VERBOSITY > 0){
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

     if(VERBOSITY > 1){
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

   if(VERBOSITY > 1) {
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

  double n_risk=0, g_risk=0, observed=0, expected=0, V=0;
  double temp1, temp2, n_events;

  vec y_time = y_node.unsafe_col(0);
  vec y_status = y_node.unsafe_col(1);

  bool break_loop = false;

  uword i = y_node.n_rows-1;

  if(VERBOSITY > 1){
   Rcout << "N obs, right node: " << sum(g_node % w_node);
   Rcout << std::endl;
  }

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
   // but turns out its faster to multiply by 0 than
   // it is to check whether n_events is > 0

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

 void Tree::grow(){

  // create inbag views of x, y, and w,
  sample_rows();

  uword n_inbag = x_inbag.n_rows;

  // all inbag observations are in the first node
  rows_node = regspace<uvec>(0, n_inbag-1);

  // assign all inbag observations to node 0
  node_assignments.zeros(n_inbag);

  // coordinate the order that nodes are grown.
  uvec nodes_open(1, fill::zeros);
  uvec nodes_queued;

  // iterate through nodes to be grown
  uvec::iterator node;


  for(node = nodes_open.begin(); node != nodes_open.end(); ++node){

   if(nodes_open[0] != 0){

    // identify which rows are in the current node.
    rows_node = find(node_assignments == *node);

   }

   y_node = y_inbag.rows(rows_node);
   w_node = w_inbag(rows_node);


   sample_cols();

   x_node = x_inbag(rows_node, cols_node);

   print_mat(x_node, "x_node", 5, 5);

   vec beta = coxph_fit(x_node,
                        y_node,
                        w_node,
                        lincomb_scale,
                        lincomb_ties_method,
                        lincomb_eps,
                        lincomb_iter_max,
                        'A');

   vec lincomb = x_node * beta;

   uvec lincomb_sort = sort_index(lincomb);

   uvec cutpoint_indices = find_cutpoints(lincomb, lincomb_sort);

   uvec cuts = linspace<uvec>(cutpoint_indices.front(),
                              cutpoint_indices.back(),
                              split_max_cuts);

   uword start = 0;

   g_node.zeros(lincomb.size());

   uvec::iterator it, it_best;

   double stat, stat_best = 0;

   for(it = cuts.begin(); it != cuts.end(); ++it){

    // flip node assignments from left to right, up to the next cutpoint
    g_node.elem(lincomb_sort.subvec(start, *it)).fill(1);

    // compute split statistics with this cut-point
    double stat = score_logrank();

    // update leaderboard
    if(stat > stat_best) {
     stat_best = stat;
     it_best = it;
    }

    start = *it;

   }




   // update tree parameters
   coef_values.push_back(beta);
   coef_indices.push_back(cols_node);

  }




 } // Tree::grow

 } // namespace aorsf

