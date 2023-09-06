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

 double TreeSurvival::compute_max_leaves(){

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


 } // namespace aorsf

