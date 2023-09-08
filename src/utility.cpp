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


 }
