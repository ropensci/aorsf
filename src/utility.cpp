/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "Utility.h"
#include "globals.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 void print_mat(arma::mat& x,
                std::string label,
                arma::uword max_cols,
                arma::uword max_rows){

  uword nrow_print = max_rows-1;
  uword ncol_print = max_cols-1;

  if(x.n_rows < max_rows) nrow_print = x.n_rows-1;
  if(x.n_cols < max_cols) ncol_print = x.n_cols-1;


  Rcout << "   ---- view of " << label << " ---- " << std::endl << std::endl;
  Rcout << x.submat(0, 0, nrow_print, ncol_print);
  Rcout << std::endl << std::endl;

 }

 void print_vec(arma::vec& x,
                std::string label,
                arma::uword max_elem){

  uword n_print = max_elem-1;

  if(x.size() < n_print) n_print = x.size()-1;

  Rcout << "   ---- view of " << label << " ---- " << std::endl << std::endl;
  Rcout << x.subvec(0, n_print);
  Rcout << std::endl << std::endl;

 }

 void print_uvec(arma::uvec& x,
                 std::string label,
                 arma::uword max_elem){

  vec x_vec = conv_to<vec>::from(x);

  print_vec(x_vec, label, max_elem);

 }

 void equalSplit(std::vector<uint>& result, uint start, uint end, uint num_parts) {

  result.reserve(num_parts + 1);

  // Return range if only 1 part
  if (num_parts == 1) {
   result.push_back(start);
   result.push_back(end + 1);
   return;
  }

  // Return vector from start to end+1 if more parts than elements
  if (num_parts > end - start + 1) {
   for (uint i = start; i <= end + 1; ++i) {
    result.push_back(i);
   }
   return;
  }

  uint length = (end - start + 1);
  uint part_length_short = length / num_parts;
  uint part_length_long = (uint) ceil(length / ((double) num_parts));
  uint cut_pos = length % num_parts;

  // Add long ranges
  for (uint i = start; i < start + cut_pos * part_length_long; i = i + part_length_long) {
   result.push_back(i);
  }

  // Add short ranges
  for (uint i = start + cut_pos * part_length_long; i <= end + 1; i = i + part_length_short) {
   result.push_back(i);
  }
 }

 vec find_unique_event_times(mat& y){

  vec result(y.n_rows);

  // times should already be sorted from low to high
  vec y_time = y.unsafe_col(0);
  vec y_status = y.unsafe_col(1);
  uword person = 0, i = 1;

  // find the first event time
  while(y_status[person] == 0){person++;}
  // assign it to first value in result
  result[0] = y_time[person];

  // find the rest (uniques only)
  for( ; person < y.n_rows; person++){

   // loop refers to i-1 in result, so i=0 was manually done above
   if(result[i-1] != y_time[person] && y_status[person] == 1){
    result[i] = y_time[person];
    i++;
   }

  }

  // resize preserves data (set_size does not)
  result.resize(i);

  return(result);

 }

 std::string uintToString(uint number) {
  return std::to_string(number);
 }

 std::string beautifyTime(uint seconds) {

  std::string result;

  // Add seconds, minutes, hours, days if larger than zero
  uint out_seconds = (uint) seconds % 60;
  result = uintToString(out_seconds) + " seconds";
  uint out_minutes = (seconds / 60) % 60;
  if (seconds / 60 == 0) {
   return result;
  } else if (out_minutes == 1) {
   result = "1 minute, " + result;
  } else {
   result = uintToString(out_minutes) + " minutes, " + result;
  }
  uint out_hours = (seconds / 3600) % 24;
  if (seconds / 3600 == 0) {
   return result;
  } else if (out_hours == 1) {
   result = "1 hour, " + result;
  } else {
   result = uintToString(out_hours) + " hours, " + result;
  }
  uint out_days = (seconds / 86400);
  if (out_days == 0) {
   return result;
  } else if (out_days == 1) {
   result = "1 day, " + result;
  } else {
   result = uintToString(out_days) + " days, " + result;
  }
  return result;
 }

 double compute_cstat(arma::mat& y,
                      arma::vec& w,
                      arma::vec& p,
                      bool pred_is_risklike){

  vec y_time   = y.unsafe_col(0);
  vec y_status = y.unsafe_col(1);

  uvec events = find(y_status == 1);

  // protection from case where there are no comparables.
  double total=0, concordant=0;

  for (uvec::iterator event = events.begin(); event < events.end(); ++event) {

   for(uword i = *event; i < y.n_rows; ++i){

    // tied events not counted
    if (y_time[i] > y_time[*event] || y_status[i] == 0) {

     total += ( (w[i] + w[*event]) / 2 );

     if (p[i] < p[*event]){

      concordant += ( (w[i]+w[*event]) / 2 );

     } else if (p[i] == p[*event]){

      concordant += ( (w[i]+w[*event]) / 4 );

     }

    }

   }

  }

  // it's possible there won't be any valid comparisons, so:
  if(total == 0) return(0.5);

  // code above assumes higher predictions mean more risk,
  if(pred_is_risklike) return(concordant / total);
  // if that isn't true (e.g., a survival prediction):
  return(1 - (concordant / total));

 }


 double compute_cstat(arma::mat& y,
                      arma::vec& w,
                      arma::uvec& g,
                      bool pred_is_risklike){

  // note: g must have only values of 0 and 1 to use this.
  // note: this is a little different in its approach than
  //       the analogous function for vec g. The reason it
  //       is different is that I've benchmarked these across
  //       big and small data and this version works best for
  //       uvec g while the analogous approach works best for
  //       vec g.

  vec y_time   = y.unsafe_col(0);
  vec y_status = y.unsafe_col(1);

  double total=0, concordant=0;

  for (uword i = 0; i < y.n_rows; ++i) {

   if(y_status[i] == 1){

    bool g_0 = g[i] == 0;

    for(uword j = i; j < y.n_rows; ++j){
     // ties not counted
     if (y_time[j] > y_time[i] || y_status[j] == 0) {

      total += ( (w[i]+w[j]) / 2 );

      // time_i < time_j, and person i had an event,
      // => if risk_i > risk_j we are concordant.
      // if risk_i is 0, risk_j cannot be less than risk_i
      // => best case scenario is a tie, i.e., g[j] == 0
      // => if g[j] is 0, we want to add 1/2 to concordant
      // => if g[j] is 1, we want to do nothing
      // => subtract 1 from g, multiply by -1, divide by 2
      if(g_0){

       if(g[j] == 0) concordant += ( (w[i]+w[j]) / 4 );

      } else if (g[j] == 1){

       // if risk_i is 1 and risk_j is 1, a tie
       concordant += ( (w[i]+w[j]) / 4 );

      } else {

       // if risk_i is 1 and risk_j is 0, concordance
       concordant += ( (w[i]+w[j]) / 2 );

      }

     }

    }

   }

  }

  // it's possible there won't be any valid comparisons, so:
  if(total == 0) return(0.5);
  // code above assumes higher predictions mean more risk,
  if(pred_is_risklike) return(concordant / total);
  // if that isn't true (e.g., a survival prediction):
  return(1 - (concordant / total));

 }

 }
