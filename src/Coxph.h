/*-----------------------------------------------------------------------------
 This file is part of aorsf.
Author: Byron C Jaeger
aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <armadillo>
#include "globals.h"

#ifndef COXPH_H
#define COXPH_H

 extern double temp1, temp2;
 extern arma::mat x_transforms;
 extern arma::uword n_vars, person, iter, i, j, k;

 namespace aorsf {

 // ----------------------------------------------------------------------------
 // ---------------------------- scaling functions -----------------------------
 // ----------------------------------------------------------------------------

 // scale observations in predictor matrix
 //
 // @description this scales inputs in the same way as
 //   the survival::coxph() function. The main reasons we do this
 //   are to avoid exponential overflow and to prevent the scale
 //   of inputs from impacting the estimated beta coefficients.
 //   E.g., you can try multiplying numeric inputs by 100 prior
 //   to calling orsf() with orsf_control_fast(do_scale = FALSE)
 //   and you will see that you get back a different forest.
 //
 // @param x_node matrix of predictors
 // @param w_node replication weights
 // @param x_transforms matrix used to store the means and scales
 //
 // @return modified x_node and x_transform filled with values
 //
 void x_node_scale(){

  // set aside memory for outputs
  // first column holds the mean values
  // second column holds the scale values

  x_transforms.zeros(n_vars, 2);
  vec means  = x_transforms.unsafe_col(0);   // Reference to column 1
  vec scales = x_transforms.unsafe_col(1);   // Reference to column 2

  w_node_sum = arma::sum(w_node);

  for(i = 0; i < n_vars; i++) {

   means.at(i) = arma::sum( w_node % x_node.col(i) ) / w_node_sum;

   x_node.col(i) -= means.at(i);

   scales.at(i) = arma::sum(w_node % abs(x_node.col(i)));

   if(scales.at(i) > 0)
    scales.at(i) = w_node_sum / scales.at(i);
   else
    scales.at(i) = 1.0; // rare case of constant covariate;

   x_node.col(i) *= scales.at(i);

  }

 }

 // same as above function, but just the means
 // (currently not used)
 void x_node_means(){

  x_transforms.zeros(n_vars, 1);
  w_node_sum = arma::sum(w_node);

  for(i = 0; i < n_vars; i++) {

   x_transforms.at(i, 0) = arma::sum( w_node % x_node.col(i) ) / w_node_sum;

  }

 }

 // ----------------------------------------------------------------------------
 // ---------------------------- cholesky functions ----------------------------
 // ----------------------------------------------------------------------------

 // cholesky decomposition
 //
 // @description this function is copied from the survival package and
 //    translated into arma.
 //
 // @param vmat matrix with covariance estimates
 // @param n_vars the number of predictors used in the current node
 //
 // prepares vmat for cholesky_solve()

 void cholesky(){

  double eps_chol = 0;
  double toler = 1e-8;
  double pivot;

  for(i = 0; i < n_vars; i++){

   if(vmat.at(i,i) > eps_chol) eps_chol = vmat.at(i,i);

   // copy upper right values to bottom left
   for(j = (i+1); j<n_vars; j++){
    vmat.at(j,i) = vmat.at(i,j);
   }
  }

  if (eps_chol == 0)
   eps_chol = toler; // no positive diagonals!
  else
   eps_chol = eps_chol * toler;

  for (i = 0; i < n_vars; i++) {

   pivot = vmat.at(i, i);

   if (pivot < Rcpp::R_PosInf && pivot > eps_chol) {

    for(j = (i+1); j < n_vars; j++){

     temp1 = vmat.at(j,i) / pivot;
     vmat.at(j,i) = temp1;
     vmat.at(j,j) -= temp1*temp1*pivot;

     for(k = (j+1); k < n_vars; k++){

      vmat.at(k, j) -= temp1 * vmat.at(k, i);

     }

    }

   } else {

    vmat.at(i, i) = 0;

   }

  }

 }

 // solve cholesky decomposition
 //
 // @description this function is copied from the survival package and
 //   translated into arma. Prepares u, the vector used to update beta.
 //
 // @param vmat matrix with covariance estimates
 // @param n_vars the number of predictors used in the current node
 //
 //
 void cholesky_solve(){

  for (i = 0; i < n_vars; i++) {

   temp1 = u[i];

   for (j = 0; j < i; j++){

    temp1 -= u[j] * vmat.at(i, j);
    u[i] = temp1;

   }

  }


  for (i = n_vars; i >= 1; i--){

   if (vmat.at(i-1, i-1) == 0){

    u[i-1] = 0;

   } else {

    temp1 = u[i-1] / vmat.at(i-1, i-1);

    for (j = i; j < n_vars; j++){
     temp1 -= u[j] * vmat.at(j, i-1);
    }

    u[i-1] = temp1;

   }

  }

 }

 // invert the cholesky in the lower triangle
 //
 // @description this function is copied from the survival package and
 //   translated into arma. Inverts vmat
 //
 // @param vmat matrix with covariance estimates
 // @param n_vars the number of predictors used in the current node
 //

 void cholesky_invert(){

  for (i=0; i<n_vars; i++){

   if (vmat.at(i,i) >0) {

    // take full advantage of the cholesky's diagonal of 1's
    vmat.at(i,i) = 1.0 / vmat.at(i,i);

    for (j=(i+1); j<n_vars; j++) {

     vmat.at(j, i) = -vmat.at(j, i);

     for (k=0; k<i; k++){
      vmat.at(j, k) += vmat.at(j, i) * vmat.at(i, k);
     }

    }

   }

  }

  /*
   ** lower triangle now contains inverse of cholesky
   ** calculate F'DF (inverse of cholesky decomp process) to get inverse
   **   of original vmat
   */
  for (i=0; i<n_vars; i++) {

   if (vmat.at(i, i) == 0) {

    for (j=0; j<i; j++) vmat.at(i, j) = 0;
    for (j=i; j<n_vars; j++) vmat.at(j, i) = 0;

   } else {

    for (j=(i+1); j<n_vars; j++) {

     temp1 = vmat.at(j, i) * vmat.at(j, j);

     if (j!=i) vmat.at(i, j) = temp1;

     for (k=i; k<j; k++){
      vmat.at(i, k) += temp1*vmat.at(j, k);
     }

    }

   }

  }

 }


 // ----------------------------------------------------------------------------
 // ------------------- Newton Raphson algo for Cox PH model -------------------
 // ----------------------------------------------------------------------------


 // Iterate (past 1) the newton raphson procedure
 //
 // @description once the first iteration is done, this procedure
 //   is used to continue updating beta
 //
 // @param beta the vector of beta coefficients
 // @param x_node the predictor matrix for the current node
 // @param y_node the outcome matrix for the current node
 // @param w_node the weight vector for the current node
 // @param n_vars the number of predictors in the current node
 // @param cmat sums of squares from x for non-events
 // @param cmat2 sums of squares from x for events
 // @param vmat matrix of covariance estimates
 // @param a sums of x for non-events
 // @param a2 sums of x for events
 // @param u vector of updates for beta
 // @param cph_method how to handle ties
 //
 // @returns the partial log likelihood based on beta
 //


 double newtraph_cph_iter(const arma::vec& beta){

  denom = 0;

  loglik = 0;

  n_risk = 0;

  person = x_node.n_rows - 1;

  u.fill(0);
  a.fill(0);
  a2.fill(0);
  vmat.fill(0);
  cmat.fill(0);
  cmat2.fill(0);

  // this loop has a strange break condition to accomodate
  // the restriction that a uvec (uword) cannot be < 0

  break_loop = false;

  XB = x_node * beta;
  Risk = exp(XB) % w_node;

  for( ; ; ){

   temp2 = y_node.at(person, 0); // time of event for current person
   n_events  = 0 ; // number of deaths at this time point
   weight_events = 0 ; // sum of w_node for the deaths
   denom_events = 0 ; // sum of weighted risks for the deaths

   // walk through this set of tied times
   while(y_node.at(person, 0) == temp2){

    n_risk++;

    xb = XB.at(person);
    risk = Risk.at(person);

    // xb = 0;
    //
    // for(i = 0; i < n_vars; i++){
    //   xb += beta.at(i) * x_node.at(person, i);
    // }

    w_node_person = w_node.at(person);

    // risk = exp(xb) * w_node_person;

    if (y_node.at(person, 1) == 0) {

     denom += risk;

     /* a contains weighted sums of x, cmat sums of squares */

     for (i=0; i<n_vars; i++) {

      temp1 = risk * x_node.at(person, i);

      a[i] += temp1;

      for (j = 0; j <= i; j++){
       cmat.at(j, i) += temp1 * x_node.at(person, j);
      }

     }

    } else {

     n_events++;

     weight_events += w_node_person;
     denom_events += risk;
     loglik += w_node_person * xb;

     for (i=0; i<n_vars; i++) {

      u[i]  += w_node_person * x_node.at(person, i);
      a2[i] += risk * x_node.at(person, i);

      for (j=0; j<=i; j++){
       cmat2.at(j, i) += risk * x_node.at(person, i) * x_node.at(person, j);
      }

     }

    }

    if(person == 0){
     break_loop = true;
     break;
    }

    person--;

   }

   // we need to add to the main terms
   if (n_events > 0) {

    if (cph_method == 0 || n_events == 1) { // Breslow

     denom  += denom_events;
     loglik -= weight_events * log(denom);

     for (i=0; i<n_vars; i++) {

      a[i]  += a2[i];
      temp1  = a[i] / denom;  // mean
      u[i]  -=  weight_events * temp1;

      for (j=0; j<=i; j++) {
       cmat.at(j, i) += cmat2.at(j, i);
       vmat.at(j, i) += weight_events * (cmat.at(j, i) - temp1 * a[j]) / denom;
      }

     }

    } else {
     /* Efron
      **  If there are 3 deaths we have 3 terms: in the first the
      **  three deaths are all in, in the second they are 2/3
      **  in the sums, and in the last 1/3 in the sum.  Let k go
      **  1 to n_events: we sequentially add a2/n_events and cmat2/n_events
      **  and efron_wt/n_events to the totals.
      */
     weight_avg = weight_events/n_events;

     for (k=0; k<n_events; k++) {

      denom  += denom_events / n_events;
      loglik -= weight_avg * log(denom);

      for (i=0; i<n_vars; i++) {

       a[i] += a2[i] / n_events;
       temp1 = a[i]  / denom;
       u[i] -= weight_avg * temp1;

       for (j=0; j<=i; j++) {
        cmat.at(j, i) += cmat2.at(j, i) / n_events;
        vmat.at(j, i) += weight_avg * (cmat.at(j, i) - temp1 * a[j]) / denom;
       }

      }

     }

    }

    a2.fill(0);
    cmat2.fill(0);

   }

   if(break_loop) break;

  }

  return(loglik);

 }

 // Do first iteration of the newton raphson procedure
 // same as above, optimized for starting value of beta = 0

 double newtraph_cph_init(){

  denom = 0;
  loglik = 0;
  n_risk = 0;

  person = n_rows - 1;

  u.fill(0);
  a.fill(0);
  a2.fill(0);
  vmat.fill(0);
  cmat.fill(0);
  cmat2.fill(0);

  // this loop has a strange break condition to accomodate
  // the restriction that a uvec (uword) cannot be < 0

  break_loop = false;

  // xb = 0.0;

  for( ; ; ){

   temp2 = y_node.at(person, 0); // time of event for current person
   n_events  = 0 ; // number of deaths at this time point
   weight_events = 0 ; // sum of w_node for the deaths
   denom_events = 0 ; // sum of weighted risks for the deaths

   // walk through this set of tied times
   while(y_node.at(person, 0) == temp2){

    n_risk++;

    risk = w_node.at(person);

    if (y_node.at(person, 1) == 0) {

     denom += risk;

     /* a contains weighted sums of x, cmat sums of squares */

     for (i=0; i<n_vars; i++) {

      temp1 = risk * x_node.at(person, i);

      a[i] += temp1;

      for (j = 0; j <= i; j++){
       cmat.at(j, i) += temp1 * x_node.at(person, j);
      }

     }

    } else {

     n_events++;

     denom_events += risk;
     // loglik += risk * xb;

     for (i=0; i<n_vars; i++) {

      temp1 = risk * x_node.at(person, i);

      u[i]  += temp1;
      a2[i] += temp1;

      for (j=0; j<=i; j++){
       cmat2.at(j, i) += temp1 * x_node.at(person, j);
      }

     }

    }

    if(person == 0){
     break_loop = true;
     break;
    }

    person--;

   }

   // we need to add to the main terms
   if (n_events > 0) {

    if (cph_method == 0 || n_events == 1) { // Breslow

     denom  += denom_events;
     loglik -= denom_events * log(denom);

     for (i=0; i<n_vars; i++) {

      a[i]  += a2[i];
      temp1  = a[i] / denom;  // mean
      u[i]  -=  denom_events * temp1;

      for (j=0; j<=i; j++) {
       cmat.at(j, i) += cmat2.at(j, i);
       vmat.at(j, i) += denom_events * (cmat.at(j, i) - temp1 * a[j]) / denom;
      }

     }

    } else {
     /* Efron
      **  If there are 3 deaths we have 3 terms: in the first the
      **  three deaths are all in, in the second they are 2/3
      **  in the sums, and in the last 1/3 in the sum.  Let k go
      **  1 to n_events: we sequentially add a2/n_events and cmat2/n_events
      **  and efron_wt/n_events to the totals.
      */
     weight_avg = denom_events/n_events;

     for (k = 0; k < n_events; k++) {

      denom  += denom_events / n_events;
      loglik -= weight_avg * log(denom);

      for (i = 0; i < n_vars; i++) {

       a[i] += a2[i] / n_events;
       temp1 = a[i]  / denom;
       u[i] -= weight_avg * temp1;

       for (j=0; j<=i; j++) {
        cmat.at(j, i) += cmat2.at(j, i) / n_events;
        vmat.at(j, i) += weight_avg * (cmat.at(j, i) - temp1 * a[j]) / denom;
       }

      }

     }

    }

    a2.fill(0);
    cmat2.fill(0);

   }

   if(break_loop) break;

  }

  return(loglik);

 }

 // run the newton raphson procedure
 //
 // @description identify a linear combination of predictors.
 //   This function is copied from the survival package and
 //   translated into arma with light modifications for efficiency.
 //   The procedure works with the partial likelihood function
 //   of the Cox model. All inputs are described above
 //   in newtraph_cph_iter()
 arma::vec newtraph_cph(arma::mat& x_node,
                        arma::mat& y_node,
                        arma::vec& w_node,
                        arma::vec& vi_pval_numer,
                        arma::vec& vi_pval_denom,
                        bool cph_do_scale,
                        int cph_method,
                        double cph_eps,
                        int cph_iter_max,
                        VariableImportance variable_importance){


  arma::uword n_vars = x_node.n_cols;
  arma::uword n_rows = x_node.n_rows;
  arma::uvec cols_node = arma::regspace<arma::uvec>(0, n_vars - 1);
  arma::vec beta_current(n_vars, arma::fill::zeros);
  arma::vec beta_new(n_vars, arma::fill::zeros);

  double temp1, temp2, halving, stat_best, stat_current;

  // these are filled with initial values later
  arma::vec XB(n_rows);
  arma::vec Risk(n_rows);
  arma::vec u(n_vars);
  arma::vec a(n_vars);
  arma::vec a2(n_vars);
  arma::mat vmat(n_vars, n_vars);
  arma::mat cmat(n_vars, n_vars);
  arma::mat cmat2(n_vars, n_vars);

  x_node_scale();

  halving = 0;

  // do the initial iteration
  stat_best = newtraph_cph_init();

  // update beta_current
  cholesky();
  cholesky_solve();
  beta_new = beta_current + u;

  if(cph_iter_max > 1 && stat_best < R_PosInf){

   for(iter = 1; iter < cph_iter_max; iter++){

    // if(verbose > 0){
    //
    //  Rcout << "--------- Newt-Raph algo; iter " << iter;
    //  Rcout << " ---------"  << std::endl;
    //  Rcout << "beta: "      << beta_new.t();
    //  Rcout << "loglik:    " << stat_best;
    //  Rcout                  << std::endl;
    //  Rcout << "------------------------------------------";
    //  Rcout << std::endl << std::endl << std::endl;
    //
    // }

    // do the next iteration
    stat_current = newtraph_cph_iter(beta_new);

    cholesky();

    // don't go trying to fix this, just use the last
    // set of valid coefficients
    if(std::isinf(stat_current)) break;

    // check for convergence
    // break the loop if the new ll is ~ same as old best ll
    if(fabs(1 - stat_best / stat_current) < cph_eps){
     break;
    }

    if(stat_current < stat_best){ // it's not converging!

     halving++; // get more aggressive when it doesn't work

     // reduce the magnitude by which beta_new modifies beta_current
     for (i = 0; i < n_vars; i++){
      beta_new[i] = (beta_new[i]+halving*beta_current[i]) / (halving+1.0);
     }

     // yeah its not technically the best but I need to do this for
     // more reasonable output when verbose = true; I should remove
     // this line when verbosity is taken out.
     stat_best = stat_current;

    } else { // it's converging!

     halving = 0;
     stat_best = stat_current;

     cholesky_solve();

     for (i = 0; i < n_vars; i++) {

      beta_current[i] = beta_new[i];
      beta_new[i] = beta_new[i] +  u[i];

     }

    }

   }

  }

  // invert vmat
  cholesky_invert();

  for (i=0; i < n_vars; i++) {

   beta_current[i] = beta_new[i];

   if(std::isinf(beta_current[i]) || std::isnan(beta_current[i])){
    beta_current[i] = 0;
   }

   if(std::isinf(vmat.at(i, i)) || std::isnan(vmat.at(i, i))){
    vmat.at(i, i) = 1.0;
   }

   // if(verbose > 0) Rcout << "scaled beta: " << beta_current[i] << "; ";

   if(cph_do_scale){
    beta_current.at(i) *= x_transforms.at(i, 1);
    vmat.at(i, i) *= x_transforms.at(i, 1) * x_transforms.at(i, 1);
   }

   // if(verbose > 0) Rcout << "un-scaled beta: " << beta_current[i] << std::endl;

   if(variable_importance == VI_ANOVA){

    if(beta_current.at(i) != 0){

     temp1 = R::pchisq(pow(beta_current[i], 2) / vmat.at(i, i),
                       1, false, false);

     if(temp1 < 0.01) vi_pval_numer[cols_node[i]]++;

    }

    vi_pval_denom[cols_node[i]]++;

   }

  }

  // if(verbose > 1) Rcout << std::endl;

  return(beta_current);

 }

 }

#endif /* COXPH_H */

