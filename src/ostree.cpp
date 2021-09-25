#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// ----------------------------------------------------------------------------
// ---------------------------- global parameters -----------------------------
// ----------------------------------------------------------------------------

// special note: dont change these doubles to uword,
//               even though some of them could be uwords;
//               operations involving uwords and doubles are not
//               straightforward and may break the routine.
// also: double + uword is slower than double + double.

double
 weight_avg,
 weight_events,
 denom_events,
 denom,
 n_events,
 n_at_risk,
 person_time,
 w_node_person,
 x_beta,
 risk,
 loglik,
 temp;

// armadillo unsigned integers
arma::uword
 i,
 j,
 k,
 iter,
 person,
 n_vars;

// a delayed break statement
bool break_loop;

// armadillo vectors (doubles)
arma::vec
 beta_current,
 beta_new,
 w_node,
 u,
 a,
 a2,
 XB,
 Risk;

// armadillo matrices (doubles)
arma::mat
 x_node,
 y_node,
 imat,
 cmat,
 cmat2;

// [[Rcpp::export]]
arma::mat x_scale_wtd(){

 // set aside memory for outputs
 // first column holds the mean values
 // second column holds the scale values

 arma::mat out(n_vars, 2);
 arma::vec means = out.unsafe_col(0);   // Reference to column 1
 arma::vec scales = out.unsafe_col(1);  // Reference to column 2

 double w_node_sum = arma::sum(w_node);

 for(i = 0; i < n_vars; i++) {

  arma::vec x_i = x_node.unsafe_col(i);

  means.at(i) = arma::sum( w_node % x_i ) / w_node_sum;

  x_i -= means.at(i);

  scales.at(i) = arma::sum(w_node % arma::abs(x_i));

  if(scales(i) > 0)
   scales.at(i) = w_node_sum / scales.at(i);
  else
   scales.at(i) = 1.0; // rare case of constant covariate;

  x_i *= scales.at(i);

 }


 return(out);

}

// ----------------------------------------------------------------------------
// ---------------------------- cholesky functions ----------------------------
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
void cholesky(){

 double eps_chol = 0;
 double toler = 1e-8;
 double pivot;

 for(i = 0; i < n_vars; i++){

  if(imat.at(i,i) > eps_chol) eps_chol = imat.at(i,i);

  // copy upper right values to bottom left
  for(j = (i+1); j<n_vars; j++){
   imat.at(j,i) = imat.at(i,j);
  }
 }

 if (eps_chol == 0)
  eps_chol = toler; // no positive diagonals!
 else
  eps_chol = eps_chol * toler;

 for (i = 0; i < n_vars; i++) {

  pivot = imat.at(i, i);

  if (pivot < R_PosInf && pivot > eps_chol) {

   for(j = (i+1); j < n_vars; j++){

    temp = imat.at(j,i) / pivot;
    imat.at(j,i) = temp;
    imat.at(j,j) -= temp*temp*pivot;

    for(k = (j+1); k < n_vars; k++){

     imat.at(k, j) -= temp * imat.at(k, i);

    }

   }

  } else {

   imat.at(i, i) = 0;

  }

 }

}

// [[Rcpp::export]]
void cholesky_solve(){

 for (i = 0; i < n_vars; i++) {

  temp = u[i];

  for (j = 0; j < i; j++){

   temp -= u[j] * imat.at(i, j);
   u[i] = temp;

  }

 }


 for (i = n_vars; i >= 1; i--){

  if (imat.at(i-1, i-1) == 0){

   u[i-1] = 0;

  } else {

   temp = u[i-1] / imat.at(i-1, i-1);

   for (j = i; j < n_vars; j++){
    temp -= u[j] * imat.at(j, i-1);
   }

   u[i-1] = temp;

  }

 }

}

// [[Rcpp::export]]
void cholesky_invert(){

 /*
  ** invert the cholesky in the lower triangle
  **   take full advantage of the cholesky's diagonal of 1's
  */
 for (i=0; i<n_vars; i++){

  if (imat.at(i,i) >0) {

   imat.at(i,i) = 1.0 / imat.at(i,i);

   for (j=(i+1); j<n_vars; j++) {

    imat.at(j, i) = -imat.at(j, i);

    for (k=0; k<i; k++){
     imat.at(j, k) += imat.at(j, i) * imat.at(i, k);
    }

   }

  }

 }

 /*
  ** lower triangle now contains inverse of cholesky
  ** calculate F'DF (inverse of cholesky decomp process) to get inverse
  **   of original imat
  */
 for (i=0; i<n_vars; i++) {

  if (imat.at(i, i) == 0) {

   for (j=0; j<i; j++) imat.at(i, j) = 0;
   for (j=i; j<n_vars; j++) imat.at(j, i) = 0;

  } else {

   for (j=(i+1); j<n_vars; j++) {

    temp = imat.at(j, i) * imat.at(j, j);

    if (j!=i) imat.at(i, j) = temp;

    for (k=i; k<j; k++){
     imat.at(i, k) += temp*imat.at(j, k);
    }

   }

  }

 }

}

// ----------------------------------------------------------------------------
// ------------------- Newton Raphson algo for Cox PH model -------------------
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
double newtraph_cph_iter(const arma::uword& method,
                         const arma::vec& beta){

 denom = 0;

 loglik = 0;

 n_at_risk = 0;

 person = x_node.n_rows - 1;

 u.fill(0);
 a.fill(0);
 a2.fill(0);
 imat.fill(0);
 cmat.fill(0);
 cmat2.fill(0);

 // this loop has a strange break condition to accomodate
 // the restriction that a uvec (uword) cannot be < 0

 break_loop = false;

 XB(arma::span(0, person)) = x_node * beta;
 Risk(arma::span(0, person)) = arma::exp(XB) % w_node;


 for( ; ; ){

  person_time = y_node.at(person, 0); // time of event for current person
  n_events  = 0 ; // number of deaths at this time point
  weight_events = 0 ; // sum of w_node for the deaths
  denom_events = 0 ; // sum of weighted risks for the deaths

  // walk through this set of tied times
  while(y_node.at(person, 0) == person_time){

   n_at_risk++;

   x_beta = XB.at(person);
   risk = Risk.at(person);

   // x_beta = 0;
   //
   // for(i = 0; i < n_vars; i++){
   //   x_beta += beta_current.at(i) * x.at(person, i);
   // }

   w_node_person = w_node.at(person);

   //risk = exp(x_beta) * w_node_person;

   if (y_node.at(person, 1) == 0) {

    denom += risk;

    /* a contains weighted sums of x, cmat sums of squares */

    for (i=0; i<n_vars; i++) {

     temp = risk * x_node.at(person, i);

     a[i] += temp;

     for (j = 0; j <= i; j++){
      cmat.at(j, i) += temp * x_node.at(person, j);
     }

    }

   } else {

    n_events++;

    weight_events += w_node_person;
    denom_events += risk;
    loglik += w_node_person * x_beta;

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

   if (method == 0 || n_events == 1) { // Breslow

    denom  += denom_events;
    loglik -= weight_events * log(denom);

    for (i=0; i<n_vars; i++) {

     a[i]  += a2[i];
     temp  = a[i] / denom;  // mean
     u[i]  -=  weight_events * temp;

     for (j=0; j<=i; j++) {
      cmat.at(j, i) += cmat2.at(j, i);
      imat.at(j, i) += weight_events * (cmat.at(j, i) - temp * a[j]) / denom;
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
      temp = a[i]  / denom;
      u[i] -= weight_avg * temp;

      for (j=0; j<=i; j++) {
       cmat.at(j, i) += cmat2.at(j, i) / n_events;
       imat.at(j, i) += weight_avg * (cmat.at(j, i) - temp * a[j]) / denom;
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


// [[Rcpp::export]]
double newtraph_cph_init(const arma::uword& method){

 denom = 0;
 loglik = 0;
 n_at_risk = 0;

 person = x_node.n_rows - 1;

 u.fill(0);
 a.fill(0);
 a2.fill(0);
 imat.fill(0);
 cmat.fill(0);
 cmat2.fill(0);

 // this loop has a strange break condition to accomodate
 // the restriction that a uvec (uword) cannot be < 0

 break_loop = false;

 x_beta = 0.0;

 for( ; ; ){

  person_time = y_node.at(person, 0); // time of event for current person
  n_events  = 0 ; // number of deaths at this time point
  weight_events = 0 ; // sum of w_node for the deaths
  denom_events = 0 ; // sum of weighted risks for the deaths

  // walk through this set of tied times
  while(y_node.at(person, 0) == person_time){

   n_at_risk++;

   risk = w_node.at(person);

   if (y_node.at(person, 1) == 0) {

    denom += risk;

    /* a contains weighted sums of x, cmat sums of squares */

    for (i=0; i<n_vars; i++) {

     temp = risk * x_node.at(person, i);

     a[i] += temp;

     for (j = 0; j <= i; j++){
      cmat.at(j, i) += temp * x_node.at(person, j);
     }

    }

   } else {

    n_events++;

    denom_events += risk;
    loglik += risk * x_beta;

    for (i=0; i<n_vars; i++) {

     temp = risk * x_node.at(person, i);

     u[i]  += temp;
     a2[i] += temp;

     for (j=0; j<=i; j++){
      cmat2.at(j, i) += temp * x_node.at(person, j);
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

   if (method == 0 || n_events == 1) { // Breslow

    denom  += denom_events;
    loglik -= denom_events * log(denom);

    for (i=0; i<n_vars; i++) {

     a[i]  += a2[i];
     temp  = a[i] / denom;  // mean
     u[i]  -=  denom_events * temp;

     for (j=0; j<=i; j++) {
      cmat.at(j, i) += cmat2.at(j, i);
      imat.at(j, i) += denom_events * (cmat.at(j, i) - temp * a[j]) / denom;
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
      temp = a[i]  / denom;
      u[i] -= weight_avg * temp;

      for (j=0; j<=i; j++) {
       cmat.at(j, i) += cmat2.at(j, i) / n_events;
       imat.at(j, i) += weight_avg * (cmat.at(j, i) - temp * a[j]) / denom;
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


// [[Rcpp::export]]
arma::mat newtraph_cph(NumericMatrix& x,
                       NumericMatrix& y,
                       NumericVector& weights,
                       const arma::uword& method,
                       const double& eps,
                       const arma::uword& iter_max,
                       const bool& rescale){

 x_node = arma::mat(x.begin(), x.nrow(), x.ncol(), false);
 y_node = arma::mat(y.begin(), y.nrow(), y.ncol(), false);
 w_node = arma::vec(weights.begin(), weights.size(), false);

 // x_node = x;
 // y_node = y;
 // w_node = weights;

 n_vars = x_node.n_cols;

 arma::mat x_transforms = x_scale_wtd();

 double ll_new, ll_best, halving = 0;

 // set_size is fast but unsafe, initial values are random
 beta_current.set_size(n_vars);
 beta_new.set_size(n_vars);

 // fill with 0's, o.w. different betas every time you run this
 beta_current.fill(0);
 beta_new.fill(0);


 // these are filled with initial values later
 XB.set_size(x_node.n_rows);
 Risk.set_size(x_node.n_rows);
 u.set_size(n_vars);
 a.set_size(n_vars);
 a2.set_size(n_vars);
 imat.set_size(n_vars, n_vars);
 cmat.set_size(n_vars, n_vars);
 cmat2.set_size(n_vars, n_vars);

 // do the initial iteration
 ll_best = newtraph_cph_init(method);

 //Rcpp::Rcout << "ll_best: " << ll_best << std::endl;

 // update beta_current
 cholesky();
 cholesky_solve();
 beta_new = beta_current + u;

 if(iter_max > 1 && ll_best < R_PosInf){

  for(iter = 1; iter < iter_max; iter++){

   // do the next iteration
   ll_new = newtraph_cph_iter(method, beta_new);

   //Rcpp::Rcout << "ll_new: " << ll_new << std::endl;

   cholesky();

   // check for convergence
   if(abs(1 - ll_best / ll_new) < eps){
    break;
   }

   if(ll_new < ll_best){ // it's not converging!

    halving++; // get more aggressive when it doesn't work

    // reduce the magnitude by which beta_new modifies beta_current
    for (i = 0; i < n_vars; i++){
     beta_new[i] = (beta_new[i]+halving*beta_current[i]) / (halving+1.0);
    }


   } else { // it's converging!

    halving = 0;
    ll_best = ll_new;

    cholesky_solve();

    for (i = 0; i < n_vars; i++) {

     beta_current[i] = beta_new[i];
     beta_new[i] = beta_new[i] +  u[i];

    }

   }

  }

 }

 // invert imat and return to original scale
 cholesky_invert();

 for (i = 0; i < n_vars; i++) {
  beta_current[i] = beta_new[i];
 }

 if(rescale == true){

  for (i=0; i < n_vars; i++) {

   beta_current[i] *= x_transforms.at(i, 1);
   u[i] /= x_transforms.at(i, 1);
   imat.at(i, i) *= x_transforms.at(i, 1) * x_transforms.at(i, 1);

   for (j = 0; j < i; j++) {

    imat.at(j, i) *= x_transforms.at(i, 1) * x_transforms.at(j, 1);
    imat.at(i, j) = imat.at(j, i);

   }

  }

 }


 for(i = 0; i < n_vars; i++){

  if(std::isinf(beta_current[i])) beta_current[i] = 0;

  if(std::isinf(imat.at(i, i))) imat.at(i, i) = 1.0;

 }

 arma::vec se = arma::sqrt(imat.diag());

 arma::vec pv(n_vars);

 for(i = 0; i < n_vars; i++){
  pv[i] = R::pchisq(pow(beta_current[i]/se[i], 2), 1, false, false);
 }

 arma::mat out(n_vars, 3);

 out.col(0) = beta_current;
 out.col(1) = se;
 out.col(2) = pv;

 return(out);

}

