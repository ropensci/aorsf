/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "globals.h"
#include "Coxph.h"
#include "utility.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 void cholesky_decomp(mat& vmat){

  double eps_chol = 0;
  double toler = 1e-8;
  double pivot1, pivot2;
  uword n_vars = vmat.n_cols;
  uword i, j, k;

  for(i = 0; i < n_vars; i++){

   if(vmat.at(i,i) > eps_chol) eps_chol = vmat.at(i,i);

   // copy upper right values to bottom left
   for(j = (i+1); j<n_vars; j++){
    vmat.at(j,i) = vmat.at(i,j);
   }
  }

  eps_chol = toler;

  for (i = 0; i < n_vars; i++) {

   pivot1 = vmat.at(i, i);

   if (pivot1 < R_PosInf && pivot1 > eps_chol) {

    for(j = (i+1); j < n_vars; j++){

     pivot2 = vmat.at(j,i) / pivot1;
     vmat.at(j,i) = pivot2;
     vmat.at(j,j) -= pivot2*pivot2*pivot1;

     for(k = (j+1); k < n_vars; k++){

      vmat.at(k, j) -= pivot2 * vmat.at(k, i);

     }

    }

   } else {

    vmat.at(i, i) = 0;

   }

  }

 }


 void cholesky_solve(mat& vmat,
                     vec& u){

  uword n_vars = vmat.n_cols;
  uword i, j;
  double temp;

  for (i = 0; i < n_vars; i++) {

   temp = u[i];

   for (j = 0; j < i; j++){

    temp -= u[j] * vmat.at(i, j);
    u[i] = temp;

   }

  }


  for (i = n_vars; i >= 1; i--){

   if (vmat.at(i-1, i-1) == 0){

    u[i-1] = 0;

   } else {

    temp = u[i-1] / vmat.at(i-1, i-1);

    for (j = i; j < n_vars; j++){
     temp -= u[j] * vmat.at(j, i-1);
    }

    u[i-1] = temp;

   }

  }

 }

 void cholesky_invert(mat& vmat){

  uword n_vars = vmat.n_cols;
  uword i, j, k;
  double temp;

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

     temp = vmat.at(j, i) * vmat.at(j, j);

     if (j!=i) vmat.at(i, j) = temp;

     for (k=i; k<j; k++){
      vmat.at(i, k) += temp*vmat.at(j, k);
     }

    }

   }

  }

 }


 arma::mat coxph_fit(arma::mat& x_node,
                     arma::mat& y_node,
                     arma::vec& w_node,
                     bool do_scale,
                     int ties_method,
                     double epsilon,
                     arma::uword iter_max){

  uword
  person,
  iter,
  i,
  j,
  k;

  vec
  beta_current,
  beta_new,
  Risk,
  u,
  a,
  a2,
  means,
  scales,
  XB;

  mat
  vmat,
  cmat,
  cmat2,
  x_transforms;

  bool break_loop;

  double
   temp1,
   temp2,
   halving,
   stat_best,
   denom,
   loglik,
   xb,
   risk,
   n_risk,
   n_events,
   weight_events,
   weight_avg,
   denom_events,
   w_node_person,
   w_node_sum,
   stat_current;

  uword n_vars = x_node.n_cols;

  if(do_scale){

   x_transforms.set_size(n_vars, 2);
   x_transforms.fill(0);

   means  = x_transforms.unsafe_col(0);   // Reference to column 1
   scales = x_transforms.unsafe_col(1);   // Reference to column 2

   w_node_sum = sum(w_node);

   for(i = 0; i < n_vars; i++) {

    means.at(i) = sum( w_node % x_node.col(i) ) / w_node_sum;

    x_node.col(i) -= means.at(i);

    scales.at(i) = sum(w_node % abs(x_node.col(i)));

    if(scales(i) > 0)
     scales.at(i) = w_node_sum / scales.at(i);
    else
     scales.at(i) = 1.0; // rare case of constant covariate;

    x_node.col(i) *= scales.at(i);

   }

  }

  beta_current.zeros(n_vars);
  beta_new.zeros(n_vars);

  // these are filled with initial values later
  Risk.set_size(x_node.n_rows);
  u.set_size(n_vars);
  a.set_size(n_vars);
  a2.set_size(n_vars);
  vmat.set_size(n_vars, n_vars);
  cmat.set_size(n_vars, n_vars);
  cmat2.set_size(n_vars, n_vars);

  halving = 0;

  // do the initial iteration
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


  // the outer loop needs to be broken when a condition occurs in
  // the inner loop - set up a bool to break the outer loop
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

    if (ties_method == 0 || n_events == 1) { // Breslow

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

  stat_best = loglik;


  // update beta_current
  cholesky_decomp(vmat);
  cholesky_solve(vmat, u);
  beta_new = beta_current + u;

  // for fast cph, returned XB needs to be computed
  if(iter_max <= 1) XB = x_node * beta_new;

  // for standard cph, iterate until convergence
  if(iter_max > 1 && stat_best < R_PosInf){

   for(iter = 1; iter < iter_max; iter++){

    // if(VERBOSITY > 1){
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

    // XB = x_node * beta_new;
    // Risk = exp(XB) % w_node;


    for( ; ; ){

     temp2 = y_node.at(person, 0); // time of event for current person
     n_events  = 0 ; // number of deaths at this time point
     weight_events = 0 ; // sum of w_node for the deaths
     denom_events = 0 ; // sum of weighted risks for the deaths

     // walk through this set of tied times
     while(y_node.at(person, 0) == temp2){

      n_risk++;

      // xb = XB.at(person);
      // risk = Risk.at(person);

      xb = 0;

      for(i = 0; i < n_vars; i++){
        xb += beta_new.at(i) * x_node.at(person, i);
      }

      w_node_person = w_node.at(person);

      risk = exp(xb) * w_node_person;

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

      if (ties_method == 0 || n_events == 1) { // Breslow

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

    stat_current = loglik;

    cholesky_decomp(vmat);

    // don't go trying to fix this, just use the last
    // set of valid coefficients
    if(std::isinf(stat_current)) break;

    // check for convergence
    // break the loop if the new ll is ~ same as old best ll
    if(fabs(1 - stat_best / stat_current) < epsilon){
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

     cholesky_solve(vmat, u);

     for (i = 0; i < n_vars; i++) {

      beta_current[i] = beta_new[i];
      beta_new[i] = beta_new[i] +  u[i];

     }

    }

   }

  }


  // invert vmat
  cholesky_invert(vmat);

  vec pvalues(beta_current.size());

  for (i=0; i < n_vars; i++) {

   beta_current[i] = beta_new[i];

   if(std::isinf(beta_current[i]) || std::isnan(beta_current[i])){
    beta_current[i] = 0;
   }

   if(std::isinf(vmat.at(i, i)) || std::isnan(vmat.at(i, i))){
    vmat.at(i, i) = 1.0;
   }

   pvalues[i] = R::pchisq(
    pow(beta_current[i], 2) / vmat.at(i, i), 1, false, false
   );

   if(do_scale){
    // return beta and variance to original scales
    beta_current.at(i) *= scales[i];
    vmat.at(i, i) *= scales[i] * scales[i];
    // same for x_node
    x_node.col(i) /= scales.at(i);
    x_node.col(i) += means.at(i);
   }

  }

  return(join_horiz(beta_current, pvalues));

 }

 }


