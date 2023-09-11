/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "TreeSurvival.h"
#include "Coxph.h"
#include "NodeSplitStats.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 TreeSurvival::TreeSurvival() { }

 TreeSurvival::TreeSurvival(double leaf_min_events,
                            double split_min_events,
                            arma::vec* unique_event_times,
                            arma::vec pred_horizon){

  this->leaf_min_events = leaf_min_events;
  this->split_min_events = split_min_events;
  this->unique_event_times = unique_event_times;
  this->pred_horizon = pred_horizon;

 }

 TreeSurvival::TreeSurvival(std::vector<double>& cutpoint,
                            std::vector<arma::uword>& child_left,
                            std::vector<arma::vec>& coef_values,
                            std::vector<arma::uvec>& coef_indices,
                            std::vector<arma::vec>& leaf_pred_indx,
                            std::vector<arma::vec>& leaf_pred_prob,
                            std::vector<arma::vec>& leaf_pred_chaz,
                            std::vector<double>& leaf_summary,
                            arma::vec pred_horizon) :
 Tree(cutpoint, child_left, coef_values, coef_indices, leaf_summary),
 leaf_pred_indx(leaf_pred_indx),
 leaf_pred_prob(leaf_pred_prob),
 leaf_pred_chaz(leaf_pred_chaz),
 pred_horizon(pred_horizon){ }

 void TreeSurvival::resize_leaves(arma::uword new_size) {

  leaf_pred_indx.resize(new_size);
  leaf_pred_prob.resize(new_size);
  leaf_pred_chaz.resize(new_size);
  leaf_summary.resize(new_size);

 }

 double TreeSurvival::compute_max_leaves(){

  n_events_inbag = sum(w_inbag % y_inbag.col(1));

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

  return(max_leaves);

 }

 bool TreeSurvival::is_col_splittable(uword j){

  uvec::iterator i;

  // initialize as 0 but do not make comparisons until x_first_value
  // is formally defined at the first instance of status == 1
  double x_first_value=0;

  bool x_first_undef = true;

  for (i = rows_node.begin(); i != rows_node.end(); ++i) {

   // if event occurred for this observation
   // column is only splittable if X is non-constant among
   // observations where an event occurred.
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

 bool TreeSurvival::is_node_splittable_internal(){

  double n_risk = sum(w_node);
  double n_events = sum(y_node.col(1) % w_node);

  return(n_events >= 2*leaf_min_events &&
         n_risk   >= 2*leaf_min_obs &&
         n_events >=   split_min_events &&
         n_risk   >=   split_min_obs);

 }

 uvec TreeSurvival::find_cutpoints(){

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
        n_risk   >= leaf_min_obs ) {

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

   if(lincomb[*it] != lincomb[*(it-1)]){

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

   output = {j};
   return(output);

  }

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

  output = join_vert(output_left, output_middle, output_right);

  return(output);

 }

 double TreeSurvival::compute_split_score(){

  double result;

  switch (split_rule) {

  case SPLIT_LOGRANK: {
   result = score_logrank();
   break;
  }

  case SPLIT_CONCORD: {
   result = score_concord();
   break;
  }

  }

  return(result);

 }

 double TreeSurvival::score_logrank(){

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

 double TreeSurvival::score_concord(){


  vec y_time = y_node.unsafe_col(0);
  vec y_status = y_node.unsafe_col(1);

  uvec event_indices = find(y_status == 1);

  uvec::iterator event;

  // protection from case where there are no comparables.
  double total=0.001, concordant=0;

  for (event = event_indices.begin(); event < event_indices.end(); ++event) {

   for(uword j = *event; j < y_node.n_rows; ++j){

    if (y_time[j] > y_time[*event]) { // ties not counted

     total += w_node[j];

     if (g_node[j] < g_node[*event]){

      concordant += w_node[j];

     } else if (g_node[j] == g_node[*event]){

      concordant += (w_node[j] / 2);

     }

    }

   }

  }

  return(concordant / total);

 }

 void TreeSurvival::node_sprout(uword node_id){

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
  // (TODO: should this case even occur? consider removing)
  if(person == y_node.n_rows){

   vec temp_surv(1, arma::fill::ones);
   vec temp_chf(1, arma::fill::zeros);

   leaf_pred_indx[node_id] = leaf_data.col(0);
   leaf_pred_prob[node_id] = temp_surv;
   leaf_pred_chaz[node_id] = temp_chf;
   leaf_summary[node_id] = 0.0;

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

  leaf_pred_indx[node_id] = leaf_data.col(0);
  leaf_pred_prob[node_id] = leaf_data.col(1);
  leaf_pred_chaz[node_id] = leaf_data.col(2);
  leaf_summary[node_id] = compute_mortality(leaf_data);

 }

 double TreeSurvival::compute_mortality(arma::mat& leaf_data){

  double result = 0;
  uword i=0, j=0;

  for( ; i < (*unique_event_times).size(); i++){

   if((*unique_event_times)[i] >= leaf_data.at(j, 0) &&
      j < (leaf_data.n_rows-1)) {j++;}

   result += leaf_data.at(j, 2);

  }

  return(result);

 }

 void TreeSurvival::predict_value(arma::mat* pred_output,
                                  arma::vec* pred_denom,
                                  char pred_type,
                                  bool oobag){

  uvec pred_leaf_sort = sort_index(pred_leaf, "ascend");

  uvec::iterator it = pred_leaf_sort.begin();

  // oobag leaf prediction has zeros for inbag rows
  if(oobag){
   while(pred_leaf(*it) == 0 && it < pred_leaf_sort.end()){
    ++it;
   }
  }

  if(it == pred_leaf_sort.end()){
   if(VERBOSITY > 0){
    Rcout << "Tree was empty, no predictions were made" << std::endl;
   }
   return;
  }

  double pred_t0;

  if(pred_type == 'S' || pred_type == 'R'){
   pred_t0 = 1;
  } else {
   pred_t0 = 0;
  }

  uword i, j;

  vec leaf_times, leaf_values;

  vec temp_vec(pred_horizon.size());
  double temp_dbl;

  do {

   uword leaf_id = pred_leaf[*it];

   // copies of leaf data using same aux memory
   leaf_times = vec(leaf_pred_indx[leaf_id].begin(),
                    leaf_pred_indx[leaf_id].size(),
                    false);

   leaf_values = vec(leaf_pred_prob[leaf_id].begin(),
                     leaf_pred_prob[leaf_id].size(),
                     false);

   if(leaf_values.is_empty()) Rcpp::stop("empty leaf");

   // don't reset i in the loop.
   // (wasteful b/c leaf_times ascend)
   i = 0;

   for(j = 0; j < pred_horizon.size(); j++){

    // t is the current prediction time
    double t = pred_horizon[j];

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

 double TreeSurvival::compute_prediction_accuracy(){


  vec y_time = y_oobag.unsafe_col(0);
  vec y_status = y_oobag.unsafe_col(1);

  uvec oobag_pred_leaf = pred_leaf(rows_oobag);

  vec mortality(rows_oobag.size());

  for(uword i = 0; i < mortality.size(); ++i){
   mortality[i] = leaf_summary[oobag_pred_leaf[i]];
  }

  Rcout << mortality << std::endl;

  // uvec event_indices = find(y_status == 1);
  //
  // uvec::iterator event;
  //
  // // protection from case where there are no comparables.
  // double total=0.001, concordant=0;
  //
  // for (event = event_indices.begin(); event < event_indices.end(); ++event) {
  //
  //  for(uword j = *event; j < y_node.n_rows; ++j){
  //
  //   if (y_time[j] > y_time[*event]) { // ties not counted
  //
  //    total += w_node[j];
  //
  //    if (g_node[j] < g_node[*event]){
  //
  //     concordant += w_node[j];
  //
  //    } else if (g_node[j] == g_node[*event]){
  //
  //     concordant += (w_node[j] / 2);
  //
  //    }
  //
  //   }
  //
  //  }
  //
  // }

  // return(concordant / total);
  return(0.0);

 }


 } // namespace aorsf

