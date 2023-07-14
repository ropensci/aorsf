/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "globals.h"
#include "NodeSplitStats.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 List lrt_multi(mat& y_node,
                mat& w_node,
                vec& XB,
                uword n_split,
                double split_min_stat,
                double leaf_min_events,
                double leaf_min_obs){

  // about this function - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // this function returns a cutpoint obtaining a local maximum
  // of the log-rank test (lrt) statistic. The default value (+Inf)
  // is really for diagnostic purposes. Put another way, if the
  // return value is +Inf (an impossible value for a cutpoint),
  // that means that we didn't find any valid cut-points and
  // the node cannot be grown with the current XB.
  //
  // if there is a valid cut-point, then the main side effect
  // of this function is to modify the group vector, which
  // will be used to assign observations to the two new nodes.
  //
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool break_loop = false;

  vec
   group(y_node.n_rows, fill::zeros),
   cutpoints_found,
   cutpoints_used(n_split),
   lrt_statistics(n_split);

  double
   stat_best = 0, // initialize at the lowest possible LRT stat value
   n_events = 0,
   n_risk = 0,
   g_risk = 0,
   stat_current,
   observed,
   expected,
   V,
   temp1,
   temp2;

  uword i, j, k, list_counter = 0;

  uvec
   jit_vals,
   // sort XB- we need to iterate over the sorted indices
   XB_sorted = sort_index(XB, "ascend");

  uvec::iterator
   XB_iter,
   jit,
   XB_iter_best;


  // group should be initialized as all 0s
  group.zeros(y_node.n_rows);

  // initialize at the lowest possible LRT stat value
  stat_best = 0;

  // sort XB to iterate over the sorted indices
  XB_sorted = sort_index(XB, "ascend");

  // unsafe_cols point to cols in y_node.
  vec y_time = y_node.unsafe_col(0);
  vec y_status = y_node.unsafe_col(1);

  // first determine the lowest value of XB that will
  // be a valid cut-point to split a node. A valid cut-point
  // is one that, if used, will result in at least leaf_min_obs
  // and leaf_min_events in both the left and right node.

  n_events = 0;
  n_risk = 0;

  if(VERBOSITY > 1){
   Rcout << "----- finding cut-point boundaries -----" << std::endl;
  }

  // Iterate through the sorted values of XB, in ascending order.

  for(XB_iter = XB_sorted.begin(); XB_iter < XB_sorted.end()-1; ++XB_iter){

   n_events += y_status[*XB_iter] * w_node[*XB_iter];
   n_risk += w_node[*XB_iter];

   // If we want to make the current value of XB a cut-point, we need
   // to make sure the next value of XB isn't equal to this current value.
   // Otherwise, we will have the same value of XB in both groups!

   if(VERBOSITY > 1){
    Rcout << XB(*XB_iter)     << " ---- ";
    Rcout << XB(*(XB_iter+1)) << " ---- ";
    Rcout << n_events     << " ---- ";
    Rcout << n_risk       << std::endl;
   }

   if(XB(*XB_iter) != XB(*(XB_iter+1))){

    if(VERBOSITY > 1){
     Rcout << "********* New cut-point here ********" << std::endl;
    }

    if( n_events >= leaf_min_events &&
        n_risk   >= leaf_min_obs) {

     if(VERBOSITY > 1){
      Rcout << std::endl;
      Rcout << "lower cutpoint: "         << XB(*XB_iter) << std::endl;
      Rcout << " - n_events, left node: " << n_events << std::endl;
      Rcout << " - n_risk, left node:   " << n_risk   << std::endl;
      Rcout << std::endl;
     }

     break;

    }

   }

  }

  if(VERBOSITY > 1){
   if(XB_iter >= XB_sorted.end()-1) {
    Rcout << "Could not find a valid lower cut-point" << std::endl;
   }
  }

  // set j to be the number of steps we have taken forward in XB
  j = XB_iter - XB_sorted.begin();

  // reset before finding the upper limit
  n_events=0;
  n_risk=0;

  // go from end to beginning to find upper cutpoint
  for(XB_iter = XB_sorted.end()-1; XB_iter >= XB_sorted.begin()+1; --XB_iter){

   n_events += y_status[*XB_iter] * w_node[*XB_iter];
   n_risk   += w_node[*XB_iter];
   group[*XB_iter] = 1;

   if(VERBOSITY > 1){
    Rcout << XB(*XB_iter)     << " ---- ";
    Rcout << XB(*(XB_iter-1)) << " ---- ";
    Rcout << n_events     << " ---- ";
    Rcout << n_risk       << std::endl;
   }

   if(XB(*XB_iter) != XB(*(XB_iter-1))){

    if(VERBOSITY > 1){
     Rcout << "********* New cut-point here ********" << std::endl;
    }

    if( n_events >= leaf_min_events &&
        n_risk   >= leaf_min_obs ) {

     // the upper cutpoint needs to be one step below the current
     // XB_iter value, because we use x <= cp to determine whether a
     // value x goes to the left node versus the right node. So,
     // if XB_iter currently points to 3, and the next value down is 2,
     // then we want to say the cut-point is 2 because then all
     // values <= 2 will go left, and 3 will go right. This matters
     // when 3 is the highest value in the vector.

     --XB_iter;

     if(VERBOSITY > 1){
      Rcout << std::endl;
      Rcout << "upper cutpoint: " << XB(*XB_iter) << std::endl;
      Rcout << " - n_events, right node: " << n_events    << std::endl;
      Rcout << " - n_risk, right node:   " << n_risk      << std::endl;
     }

     break;

    }

   }

  }

  // k = n steps from beginning of sorted XB to current XB_iter
  k = XB_iter + 1 - XB_sorted.begin();

  if(VERBOSITY > 1){
   Rcout << "----------------------------------------" << std::endl;
   Rcout << std::endl << std::endl;
   Rcout << "sorted XB: " << std::endl << XB(XB_sorted).t() << std::endl;
  }

  // initialize cut-point as the value of XB that XB_iter currently
  // points to, which is the upper cut-point of XB.
  XB_iter_best = XB_iter;

  // what happens if we don't have enough events or obs to split?
  // the first valid lower cut-point (at XB_sorted(k)) is > the first
  // valid upper cutpoint (current value of n_risk). Put another way,
  // k (the number of steps taken from beginning of the XB vec)
  // will be > n_rows - p, where the difference on the RHS is
  // telling us where we are after taking p steps from the end
  // of the XB vec. Returning the infinite cp is a red flag.

  if(VERBOSITY > 1){
   Rcout << "N steps from beginning to first cp: " << j << std::endl;
   Rcout << "N steps from beginning to last cp: " << k << std::endl;
   Rcout << "n potential cutpoints: " << k-j << std::endl;
  }

  if (j > k){

   if(VERBOSITY > 1) {
    Rcout << "Could not find a cut-point for this XB" << std::endl;
   }

   return(List::create(_["cutpoints"] = R_PosInf,
                       _["statistic"] = R_PosInf));

  }

  // set k to be number of steps b/t lower and upper cutpoint.
  k -= j;

  if(VERBOSITY > 1){
   Rcout << "----- initializing cutpoints -----" << std::endl;
  }

  // what happens if there are only 5 potential cut-points
  // but the value of n_split is > 5? We will just check out
  // the 5 valid cutpoints.
  if(k > n_split){
   // there is enough space to find n_split or more cutpoints
   jit_vals = linspace<uvec>(0, k, n_split);
  } else {
   // there is less than enough space to find n_split cutpoints,
   // so find as many as we can (i.e., k).
   jit_vals = linspace<uvec>(0, k, k);
  }

  cutpoints_found.resize( jit_vals.size() );

  if(j == 0) jit_vals(jit_vals.size()-1)--;

  for(k = 0; k < cutpoints_found.size(); k++){
   cutpoints_found(k) = XB(*(XB_iter_best - jit_vals(k)));
  }

  if(j == 0) jit_vals(jit_vals.size()-1)++;


  if(VERBOSITY > 1){

   Rcout << "cut-points chosen: ";

   Rcout << cutpoints_found.t();

   Rcout << "----------------------------------------" << std::endl <<
    std::endl << std::endl;

  }

  bool do_lrt = true;

  k = 0;
  j = 1;

  // begin outer loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for(jit = jit_vals.begin(); jit != jit_vals.end(); ++jit){


   for( ; j < *jit; j++){
    group(*XB_iter) = 1;
    --XB_iter;
   }

   // always do the test if we are on a boundary of jit
   if(jit == jit_vals.begin() ||
      jit == jit_vals.end()-1){

    do_lrt = true;

   } else {

    if( cutpoints_found(k) == cutpoints_found(k+1) ||
        cutpoints_found(k) == cutpoints_found(0)   ||
        *jit <= 1){

        Rcout << "SKIP" << std::endl;

        do_lrt = false;

    } else {

     while(XB(*XB_iter) == XB(*(XB_iter - 1))){

      group(*XB_iter) = 1;
      --XB_iter;
      ++j;

      if(VERBOSITY > 1){
       Rcout << "cutpoint dropped down one spot: ";
       Rcout << XB(*XB_iter) << std::endl;
      }

     }

     do_lrt = true;

    }

   }

   ++k;

   if(do_lrt){

    cutpoints_used(list_counter) = XB(*XB_iter);

    n_risk=0;
    g_risk=0;

    observed=0;
    expected=0;

    V=0;

    break_loop = false;

    i = y_node.n_rows-1;

    if(VERBOSITY > 1){
     Rcout << "sum(group==1): " << sum(group) << ";  ";
     Rcout << "sum(group==1 * w_node): " << sum(group % w_node);
     Rcout << std::endl;
     if(VERBOSITY > 1){
      Rcout << "group:" << std::endl;
      Rcout << group(XB_sorted).t() << std::endl;
     }
    }


    // begin inner loop  - - - - - - - - - - - - -  - - - - - - - - - - - - -
    for (; ;){

     temp1 = y_time[i];

     n_events = 0;

     for ( ; y_time[i] == temp1; i--) {

      n_risk += w_node[i];
      n_events += y_status[i] * w_node[i];
      g_risk += group[i] * w_node[i];
      observed += y_status[i] * group[i] * w_node[i];

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
    // end inner loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

    stat_current = pow(expected-observed, 2) / V;

    lrt_statistics(list_counter) = stat_current;

    list_counter++;

    if(VERBOSITY > 1){

     Rcout << "-------- log-rank test results --------" << std::endl;
     Rcout << "cutpoint: " << XB(*XB_iter)                  << std::endl;
     Rcout << "lrt stat: " << stat_current              << std::endl;
     Rcout << "---------------------------------------" << std::endl <<
      std::endl << std::endl;

    }

    if(stat_current > stat_best){
     XB_iter_best = XB_iter;
     stat_best = stat_current;
    }

   }
   // end outer loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

  }

  // if the log-rank test does not detect a difference at 0.05 alpha,
  // maybe it's not a good idea to split this node.

  // if(stat_best < 3.841459) return(
  //   List::create(_["cutpoints"] = R_PosInf,
  //                _["statistic"] = R_PosInf)
  // );

  if(VERBOSITY > 1){
   Rcout << "Best LRT stat: " << stat_best << std::endl;
  }

  // rewind XB_iter until it is back where it was when we got the
  // best lrt stat. While rewinding XB_iter, also reset the group
  // values so that group is as it was when we got the best
  // lrt stat.


  while(XB_iter <= XB_iter_best){
   group(*XB_iter) = 0;
   ++XB_iter;
  }

  return(List::create(_["cutpoints"] = cutpoints_used,
                      _["statistic"] = lrt_statistics,
                      _["cutpoints_found"] = cutpoints_found));

 }

 }

