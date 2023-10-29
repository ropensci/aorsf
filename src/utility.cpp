/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "utility.h"
#include "globals.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 // # nocov start
 void print_mat(arma::mat& x,
                std::string label,
                arma::uword max_cols,
                arma::uword max_rows){

  uword nrow_print = max_rows-1;
  uword ncol_print = max_cols-1;

  if(x.n_rows < max_rows) nrow_print = x.n_rows-1;
  if(x.n_cols < max_cols) ncol_print = x.n_cols-1;


  Rcout << "   -- " << label << std::endl << std::endl;
  Rcout << x.submat(0, 0, nrow_print, ncol_print);
  Rcout << std::endl << std::endl;

 }

 void print_vec(arma::vec& x,
                std::string label,
                arma::uword max_elem){

  uword n_print = max_elem-1;

  if(x.size() <= n_print) n_print = x.size()-1;
  Rcout << "   -- " << label << std::endl << std::endl;

  if(x.size() == 0){
   Rcout << "   empty vector";
  } else {
   Rcout << x.subvec(0, n_print).t();
  }

  Rcout << std::endl << std::endl;

 }

 void print_uvec(arma::uvec& x,
                 std::string label,
                 arma::uword max_elem){

  uword n_print = max_elem-1;

  if(x.size() <= n_print) n_print = x.size()-1;
  Rcout << "   -- " << label << std::endl << std::endl;

  if(x.size() == 0){
   Rcout << "   empty vector";
  } else {
   Rcout << x.subvec(0, n_print).t();
  }

  Rcout << std::endl << std::endl;

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

 // # nocov end

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


 double compute_logrank(arma::mat& y,
                        arma::vec& w,
                        arma::uvec& g){

  double n_risk=0, g_risk=0, observed=0, expected=0, V=0,
   temp1, temp2, n_events;

  vec y_time = y.unsafe_col(0);
  vec y_status = y.unsafe_col(1);

  bool break_loop = false;

  uword i = y.n_rows-1;

  // breaking condition of outer loop governed by inner loop
  for (; ;){

   temp1 = y_time(i);

   n_events = 0;

   for ( ; y_time[i] == temp1; i--) {

    n_risk += w[i];
    n_events += y_status[i] * w[i];
    g_risk += g[i] * w[i];
    observed += y_status[i] * g[i] * w[i];

    if(i == 0){
     break_loop = true;
     break;
    }

   }

   // should only do these calculations if n_events > 0,
   // but multiplying by 0 is usually faster than checking
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

 double compute_cstat_surv(arma::mat& y,
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


 double compute_cstat_surv(arma::mat& y,
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

 double compute_cstat_clsf(vec& y, vec& w, vec& p){

  uvec p_sort_index = sort_index(p);
  vec p_freqs(w.size());
  vec w_freqs(w.size());

  double p_current = p[p_sort_index[0]];
  uword freq_counter = 0;

  for(uword i = 0; i < w.size(); ++i){

   double p_new = p[p_sort_index[i]];

   if(p_new != p_current){
    p_freqs[freq_counter] = p_current;
    p_current = p_new;
    freq_counter++;
   }

   w_freqs[freq_counter] += w[p_sort_index[i]];

  }

  p_freqs[freq_counter] = p_current;
  freq_counter++;

  w_freqs.resize(freq_counter);
  p_freqs.resize(freq_counter);

  vec r = cumsum(w_freqs) - 0.5 * (w_freqs - 1);

  vec w_rank;

  interp1(p_freqs, r, p, w_rank);

  vec w_y = w % y;

  double n = sum(w_freqs);
  double n1 = sum(w_y);
  double mean_rank = dot(w_rank, w_y) / n1;
  double cstat = (mean_rank - (n1 + 1) * 0.5) / (n - n1);

  return(cstat);

 }

 double compute_cstat_clsf(vec& y, vec& w, uvec& g){

  double true_pos=0, true_neg=0, false_pos=0, false_neg=0;

  for(uword i = 0; i < g.size(); ++i){

   if(g[i] == 0 && y[i] == 0){

    true_neg += w[i];

   } else if(g[i] == 1 && y[i] == 1){

    true_pos += w[i];

   } else if(g[i] == 1) {

    false_pos += w[i];

   } else {

    false_neg += w[i];

   }

  }

  double sens = true_pos / (true_pos + false_pos);
  double spec = true_neg / (true_neg + false_neg);

  return(0.5 * (sens + spec));

 }

 double compute_gini(mat& y, vec& w, uvec& g){

  vec y_probs_0(y.n_cols);
  vec y_probs_1(y.n_cols);

  double n_0 = 0;
  double n_1 = 0;

  for(uword i = 0; i < y.n_rows; ++i){

   if(g[i] == 1){

    n_1 += w[i];
    y_probs_1 += (y.row(i) * w[i]);

   } else {

    n_0 += w[i];
    y_probs_0 += (y.row(i) * w[i]);

   }

  }

  y_probs_1 /= n_1;
  y_probs_0 /= n_0;

  y_probs_1 = join_vert(vec {1 - sum(y_probs_1)}, y_probs_1);
  y_probs_0 = join_vert(vec {1 - sum(y_probs_0)}, y_probs_0);

  double gini_1 = 1 - sum(y_probs_1 % y_probs_1);
  double gini_0 = 1 - sum(y_probs_0 % y_probs_0);

  double n_tot = n_1 + n_0;

  return(gini_1 * (n_1/n_tot) + gini_0 * (n_0/n_tot));

 }

 vec compute_pred_prob(mat& y, vec& w){

  double n_wtd = 0;
  vec pred_prob(y.n_cols, fill::zeros);

  for(uword i = 0; i < y.n_rows; ++i){
   n_wtd += w[i];
   for(uword j = 0; j < y.n_cols; ++j){
    pred_prob[j] += (y.at(i, j) * w[i]);
   }
  }

  pred_prob /= n_wtd;
  vec pred_0 = vec {1 - sum(pred_prob)};
  pred_prob = join_vert(pred_0, pred_prob);
  return(pred_prob);

 }

 mat expand_y_clsf(vec& y, uword n_class){

  mat out(y.n_rows, n_class - 1, fill::zeros);

  for(uword i = 0; i < y.n_rows; ++i){

   double yval = y[i];

   if(yval > 0){ out.at(i, yval-1) = 1; }

  }

  return(out);

 }

 arma::mat linreg_fit(arma::mat& x_node,
                      arma::mat& y_node,
                      arma::vec& w_node,
                      bool do_scale,
                      double epsilon,
                      arma::uword iter_max){

  // Add an intercept column to the design matrix
  vec intercept(x_node.n_rows, fill::ones);
  mat X = join_horiz(intercept, x_node);

  // for later steps we don't care about the intercept term b/c we don't
  // need it, but in this step including the intercept is important
  // for computing p-values of other regression coefficients.

  uword resid_df = X.n_rows - X.n_cols;

  vec beta = solve(X.t() * diagmat(w_node) * X, X.t() * (w_node % y_node));

  vec resid  = y_node - X * beta;

  double s2 = as_scalar(trans(resid) * (w_node % resid) / (resid_df));

  mat beta_cov = s2 * inv(X.t() * diagmat(w_node) * X);

  vec se = sqrt(diagvec(beta_cov));

  vec tscores = beta / se;

  // Calculate two-tailed p-values
  vec pvalues(X.n_cols);

  for (uword i = 0; i < X.n_cols; ++i) {

   double tstat = std::abs(tscores[i]);

   pvalues[i] = 2 * (1 - R::pt(tstat, resid_df, 1, 0));

  }


  mat result = join_horiz(beta, pvalues);

  return(result.rows(1, result.n_rows-1));

 }

 arma::mat logreg_fit(arma::mat& x_node,
                      arma::mat& y_node,
                      arma::vec& w_node,
                      bool do_scale,
                      double epsilon,
                      arma::uword iter_max){

  // Add an intercept column to the design matrix
  vec intercept(x_node.n_rows, fill::ones);
  const mat X = join_horiz(intercept, x_node);

  // for later steps we don't care about the intercept term b/c we don't
  // need it, but in this step including the intercept is important
  // for computing p-values of other regression coefficients.

  vec beta(X.n_cols, fill::zeros);
  mat hessian(X.n_cols, X.n_cols);

  for (uword iter = 0; iter < iter_max; ++iter) {

   vec eta = X * beta;

   // needs to be element-wise b/c valgrind
   // spots a possible memory leak otherwise
   for(uword i = 0; i < eta.size(); i++){
    eta[i] = exp(eta[i]);
   }

   vec pi = eta / (1 + eta);
   vec w = w_node % pi % (1 - pi);

   vec gradient = X.t() * ((y_node - pi) % w_node);
   hessian = -X.t() * diagmat(w) * X;

   beta -= solve(hessian, gradient);

   if (norm(gradient) < epsilon) {
    break;
   }
  }

  // Compute standard errors, z-scores, and p-values

  vec se = sqrt(diagvec(inv(-hessian)));
  vec zscores = beta / se;
  vec pvalues = 2 * (1 - normcdf(abs(zscores)));
  mat result = join_horiz(beta, pvalues);

  return(result.rows(1, result.n_rows-1));

 }

 arma::mat scale_x(arma::mat& x,
                   arma::vec& w){

  uword n_vars = x.n_cols;

  mat x_transforms(n_vars, 2, fill::zeros);

  vec means  = x_transforms.unsafe_col(0);   // Reference to column 1
  vec scales = x_transforms.unsafe_col(1);   // Reference to column 2

  double w_sum = sum(w);
  double m = w.size();

  for(uword i = 0; i < n_vars; i++) {

   means.at(i) = sum(w % x.col(i) ) / w_sum;

   // subtract the mean now so you don't have to do the subtraction
   // when computing standard deviation.
   x.col(i) -= means.at(i);

   scales.at(i) = sqrt(
    sum(w % pow(x.col(i), 2)) / ( (m-1) * w_sum / m )
   );

   if(scales(i) <= 0) scales.at(i) = 1.0; // constant covariate;

   x.col(i) /= scales.at(i);

  }

  return(x_transforms);

 }

 void unscale_x(arma::mat& x,
                arma::mat& x_transforms){

  uword n_vars = x.n_cols;

  vec means  = x_transforms.unsafe_col(0);   // Reference to column 1
  vec scales = x_transforms.unsafe_col(1);   // Reference to column 2

  for(uword i = 0; i < n_vars; i++){

   x.col(i) /= scales.at(i);
   x.col(i) += means.at(i);

  }

 }


 }
