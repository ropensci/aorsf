#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>

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
  w_node_sum,
  denom_events,
  denom,
  cph_eps,
  n_events,
  n_risk,
  n_risk_sub,
  g_risk,
  temp1,
  temp2,
  halving,
  stat_current,
  stat_best,
  w_node_person,
  x_beta,
  risk,
  loglik,
  cutpoint,
  observed,
  expected,
  V,
  leaf_min_obs,
  leaf_min_events,
  time_oobag,
  cph_pval_max;

int
  verbose,
  mtry_int;

// armadillo unsigned integers
arma::uword
  i,
  j,
  k,
  iter,
  mtry,
  mtry_temp,
  person,
  person_leaf,
  person_ref_index,
  n_vars,
  n_rows,
  n_slots,
  cph_method,
  cph_iter_max,
  n_split,
  nodes_max_guess,
  nodes_max_true,
  n_cols_to_sample,
  nn_left,
  leaf_node_counter,
  leaf_node_index_counter;

String
  node_name,
  person_leaf_name;

bool
  break_loop, // a delayed break statement
  oobag_pred,
  scale_fast;

// armadillo vectors (doubles)
arma::vec
  vec_temp,
  time_unique,
  node_assignments,
  nodes_grown,
  surv_oobag,
  denom_oobag,
  beta_current,
  beta_new,
  beta_cph,
  cutpoints,
  w_inbag,
  w_grown,
  w_node,
  group,
  u,
  a,
  a2,
  XB,
  Risk;

// armadillo unsigned integer vectors
arma::uvec
  iit_vals,
  jit_vals,
  rows_inbag,
  rows_oobag,
  rows_node,
  rows_leaf,
  rows_node_combined,
  cols_to_sample_01,
  cols_to_sample,
  cols_node,
  leaf_node_index,
  nodes_to_grow,
  nodes_to_grow_temp,
  obs_in_node,
  children_left,
  leaf_oobag;

// armadillo iterators for unsigned integer vectors
arma::uvec::iterator
  iit,
  iit_last,
  iit_best,
  jit,
  node;

// armadillo matrices (doubles)
arma::mat
  x_input,
  x_transforms,
  y_input,
  x_inbag,
  y_inbag,
  y_grown,
  x_oobag,
  y_oobag,
  x_node,
  y_node,
  node_sums,
  leaf,
  leaf_surv,
  imat,
  cmat,
  cmat2,
  betas,
  leaf_nodes;

arma::umat
  col_indices;



// [[Rcpp::export]]
String make_node_name(const arma::uword& part){

  std::ostringstream oss;
  oss << "node_" << part;
  String out = oss.str();
  return out;

}

// ----------------------------------------------------------------------------
// ---------------------------- scaling functions -----------------------------
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
void x_node_scale(){

  // set aside memory for outputs
  // first column holds the mean values
  // second column holds the scale values

  x_transforms.zeros(n_vars, 2);
  arma::vec means  = x_transforms.unsafe_col(0);   // Reference to column 1
  arma::vec scales = x_transforms.unsafe_col(1);   // Reference to column 2

  w_node_sum = arma::sum(w_node);

  for(i = 0; i < n_vars; i++) {

    means.at(i) = arma::sum( w_node % x_node.col(i) ) / w_node_sum;

    x_node.col(i) -= means.at(i);

    scales.at(i) = arma::sum(w_node % arma::abs(x_node.col(i)));

    if(scales(i) > 0)
      scales.at(i) = w_node_sum / scales.at(i);
    else
      scales.at(i) = 1.0; // rare case of constant covariate;

    x_node.col(i) *= scales.at(i);

  }

}

// [[Rcpp::export]]
void x_input_scale(){

  // set aside memory for outputs
  // first column holds the mean values
  // second column holds the scale values

  x_transforms.zeros(n_vars, 2);
  arma::vec means  = x_transforms.unsafe_col(0);   // Reference to column 1
  arma::vec scales = x_transforms.unsafe_col(1);   // Reference to column 2

  for(i = 0; i < n_vars; i++) {

    means.at(i) = arma::mean( x_input.col(i) );

    x_input.col(i) -= means.at(i);

    scales.at(i) = arma::sum(arma::abs(x_input.col(i)));

    if(scales(i) > 0)
      scales.at(i) = x_input.n_rows / scales.at(i);
    else
      scales.at(i) = 1.0; // rare case of constant covariate;

    x_input.col(i) *= scales.at(i);

  }

}


// ----------------------------------------------------------------------------
// -------------------------- leaf_surv functions -----------------------------
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
void leaf_surv_small(const arma::mat& y,
                     const arma::vec& w){

  time_unique.resize(y.n_rows);

  // use sorted y times to count the number of unique.
  // also define number at risk as the sum of the w
  n_slots = 1; // set at 1 to account for the first time

  // find the first unique event time
  person = 0;
  n_risk = 0;

  while(y.at(person, 1) == 0){
    n_risk += w.at(person);
    person++;
  }

  // now person should correspond to the first event time
  time_unique(0) = y.at(person, 0);  // see above
  temp2          = y.at(person, 0);

  for( ; person < y.n_rows; person++){

    if(temp2 != y.at(person,0) && y.at(person,1) == 1){

      time_unique.at(n_slots) = y.at(person,0);
      temp2 = y.at(person, 0);
      n_slots++;

    }

    n_risk += w.at(person);

  }

  // drop the extra zeros from time_unique
  time_unique = time_unique(arma::span(0, n_slots - 1));

  // reset for next loop
  person = 0; j = 0; temp1 = 1.0;

  // vec_temp.resize(n_slots);

  do {

    n_events   = 0;
    n_risk_sub = 0;
    temp2      = y.at(person, 0);

    while(y.at(person, 0) == temp2){

      n_risk_sub += w.at(person);

      if(y.at(person, 1) == 1){
        n_events += w.at(person);
      }

      if(person == y.n_rows-1) break;

      person++;

    }

    // only do km if a death was observed

    if(n_events > 0){

      temp1 = temp1 * (n_risk - n_events) / n_risk;

      leaf_nodes(leaf_node_counter, 0) = time_unique(j);
      leaf_nodes(leaf_node_counter, 1) = temp1;

      j++;
      leaf_node_counter++;

    }


    n_risk -= n_risk_sub;

  } while (j < n_slots);

  leaf_node_index(leaf_node_index_counter) = leaf_node_counter;
  leaf_node_index_counter++;
  // return(arma::join_horiz(time_unique, vec_temp));

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

        temp1 = imat.at(j,i) / pivot;
        imat.at(j,i) = temp1;
        imat.at(j,j) -= temp1*temp1*pivot;

        for(k = (j+1); k < n_vars; k++){

          imat.at(k, j) -= temp1 * imat.at(k, i);

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

    temp1 = u[i];

    for (j = 0; j < i; j++){

      temp1 -= u[j] * imat.at(i, j);
      u[i] = temp1;

    }

  }


  for (i = n_vars; i >= 1; i--){

    if (imat.at(i-1, i-1) == 0){

      u[i-1] = 0;

    } else {

      temp1 = u[i-1] / imat.at(i-1, i-1);

      for (j = i; j < n_vars; j++){
        temp1 -= u[j] * imat.at(j, i-1);
      }

      u[i-1] = temp1;

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

        temp1 = imat.at(j, i) * imat.at(j, j);

        if (j!=i) imat.at(i, j) = temp1;

        for (k=i; k<j; k++){
          imat.at(i, k) += temp1*imat.at(j, k);
        }

      }

    }

  }

}


// ----------------------------------------------------------------------------
// ------------------- Newton Raphson algo for Cox PH model -------------------
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
double newtraph_cph_iter(const arma::vec& beta){

  denom = 0;

  loglik = 0;

  n_risk = 0;

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

  XB = x_node * beta;
  Risk = arma::exp(XB) % w_node;

  for( ; ; ){

    temp2 = y_node.at(person, 0); // time of event for current person
    n_events  = 0 ; // number of deaths at this time point
    weight_events = 0 ; // sum of w_node for the deaths
    denom_events = 0 ; // sum of weighted risks for the deaths

    // walk through this set of tied times
    while(y_node.at(person, 0) == temp2){

      n_risk++;

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

      if (cph_method == 0 || n_events == 1) { // Breslow

        denom  += denom_events;
        loglik -= weight_events * log(denom);

        for (i=0; i<n_vars; i++) {

          a[i]  += a2[i];
          temp1  = a[i] / denom;  // mean
          u[i]  -=  weight_events * temp1;

          for (j=0; j<=i; j++) {
            cmat.at(j, i) += cmat2.at(j, i);
            imat.at(j, i) += weight_events * (cmat.at(j, i) - temp1 * a[j]) / denom;
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
              imat.at(j, i) += weight_avg * (cmat.at(j, i) - temp1 * a[j]) / denom;
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
double newtraph_cph_init(){

  denom = 0;
  loglik = 0;
  n_risk = 0;

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
        loglik += risk * x_beta;

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
            imat.at(j, i) += denom_events * (cmat.at(j, i) - temp1 * a[j]) / denom;
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
              imat.at(j, i) += weight_avg * (cmat.at(j, i) - temp1 * a[j]) / denom;
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
arma::vec newtraph_cph(){

  beta_current.zeros(n_vars);
  beta_new.zeros(n_vars);

  // these are filled with initial values later
  XB.set_size(x_node.n_rows);
  Risk.set_size(x_node.n_rows);
  u.set_size(n_vars);
  a.set_size(n_vars);
  a2.set_size(n_vars);
  imat.set_size(n_vars, n_vars);
  cmat.set_size(n_vars, n_vars);
  cmat2.set_size(n_vars, n_vars);

  halving = 0;

  // do the initial iteration
  stat_best = newtraph_cph_init();

  // update beta_current
  cholesky();
  cholesky_solve();
  beta_new = beta_current + u;

  if(cph_iter_max > 1 && stat_best < R_PosInf){

    for(iter = 1; iter < cph_iter_max; iter++){

      if(verbose > 0){

        Rcout << "--------- Newt-Raph algo; iter " << iter;
        Rcout << " ---------"  << std::endl;
        Rcout << "beta: "      << beta_new.t();
        Rcout << "loglik:    " << stat_best;
        Rcout                  << std::endl;
        Rcout << "------------------------------------------";
        Rcout << std::endl << std::endl << std::endl;

      }

      // do the next iteration
      stat_current = newtraph_cph_iter(beta_new);

      //Rcpp::Rcout << "stat_current: " << stat_current << std::endl;

      cholesky();

      // check for convergence
      // break the loop if the new ll is ~ same as old best ll
      if(abs(1 - stat_best / stat_current) < cph_eps){
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

  // invert imat
  cholesky_invert();

  for (i=0; i < n_vars; i++) {

    beta_current[i] = beta_new[i];

    if(std::isinf(beta_current[i])) beta_current[i] = 0;

    if(std::isinf(imat.at(i, i))) imat.at(i, i) = 1.0;

    if(verbose > 0) Rcout << "scaled beta: " << beta_current[i] << "; ";

    beta_current.at(i) *= x_transforms.at(i, 1);

    imat.at(i, i) *= x_transforms.at(i, 1) * x_transforms.at(i, 1);

    if(verbose > 0) Rcout << "un-scaled beta: " << beta_current[i] << std::endl;

    temp1 = R::pchisq(pow(beta_current[i], 2) / imat.at(i, i),
                      1, false, false);

    if(temp1 > cph_pval_max){
      beta_current[i] = 0;
      if(verbose > 0){
        Rcout<<"dropping coef "<<i<<" to 0; p = "<<temp1<<std::endl;
      }
    }

  }

  if(verbose > 0) Rcout << std::endl;

  return(beta_current);

}


// // ----------------------------------------------------------------------------
// // ---------------------------- node functions --------------------------------
// // ----------------------------------------------------------------------------

// [[Rcpp::export]]
double lrt_multi(){

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

  break_loop = false;

  // group should be initialized as all 0s
  group.zeros(y_node.n_rows);

  // initialize at the lowest possible LRT stat value
  stat_best = 0;

  // sort XB- we need to iterate over the sorted indices
  iit_vals = arma::sort_index(XB, "ascend");

  // unsafe columns point to specific cols in y_node.
  // this makes the code more readable and doesn't copy data
  arma::vec status = y_node.unsafe_col(1);
  arma::vec time = y_node.unsafe_col(0);

  // first determine the lowest value of XB that will
  // be a valid cut-point to split a node. A valid cut-point
  // is one that, if used, will result in at least leaf_min_obs
  // and leaf_min_events in both the left and right node.

  n_events = 0;
  n_risk = 0;

  if(verbose > 0){
    Rcout << "----- finding cut-point boundaries -----" << std::endl;
  }

  // Iterate through the sorted values of XB, in ascending order.

  for(iit = iit_vals.begin(); iit < iit_vals.end()-1; ++iit){

    n_events += status(*iit) * w_node(*iit);
    n_risk += w_node(*iit);

    // If we want to make the current value of XB a cut-point, we need
    // to make sure the next value of XB isn't equal to this current value.
    // Otherwise, we will have the same value of XB in both groups!

    // if(verbose > 1){
    //   Rcout << XB(*iit)     << " ---- ";
    //   Rcout << XB(*(iit+1)) << " ---- ";
    //   Rcout << n_events     << " ---- ";
    //   Rcout << n_risk       << std::endl;
    // }

    if(XB(*iit) != XB(*(iit+1))){

      // if(verbose > 1){
      //   Rcout << "********* New cut-point here ********" << std::endl;
      // }


      if( n_events >= leaf_min_events &&
          n_risk   >= leaf_min_obs) {

        if(verbose > 0){
          Rcout << std::endl;
          Rcout << "lower cutpoint: "         << XB(*iit) << std::endl;
          Rcout << " - n_events, left node: " << n_events << std::endl;
          Rcout << " - n_risk, left node:   " << n_risk   << std::endl;
          Rcout << std::endl;
        }

        break;

      }

    }

  }

  if(verbose > 0){
    if(iit >= iit_vals.end()-1) {
      Rcout << "Could not find a valid lower cut-point" << std::endl;
    }
  }


  j = iit - iit_vals.begin();

  // got to reset these before finding the upper limit
  n_events=0;
  n_risk=0;

  // do the first step in the loop manually since we need to
  // refer to iit+1 in all proceeding steps.

  for(iit = iit_vals.end()-1; iit >= iit_vals.begin()+1; --iit){

    n_events += status(*iit) * w_node(*iit);
    n_risk   += w_node(*iit);
    group(*iit) = 1;

    // if(verbose > 1){
    //   Rcout << XB(*iit)     << " ---- ";
    //   Rcout << XB(*(iit-1)) << " ---- ";
    //   Rcout << n_events     << " ---- ";
    //   Rcout << n_risk       << std::endl;
    // }

    if(XB(*iit) != XB(*(iit-1))){

      // if(verbose > 1){
      //   Rcout << "********* New cut-point here ********" << std::endl;
      // }

      if( n_events >= leaf_min_events &&
          n_risk   >= leaf_min_obs ) {

        // the upper cutpoint needs to be one step below the current
        // iit value, because we use x <= cp to determine whether a
        // value x goes to the left node versus the right node. So,
        // if iit currently points to 3, and the next value down is 2,
        // then we want to say the cut-point is 2 because then all
        // values <= 2 will go left, and 3 will go right. This matters
        // when 3 is the highest value in the vector.

        --iit;

        if(verbose > 0){
          Rcout << std::endl;
          Rcout << "upper cutpoint: " << XB(*iit) << std::endl;
          Rcout << " - n_events, right node: " << n_events    << std::endl;
          Rcout << " - n_risk, right node:   " << n_risk      << std::endl;
        }

        break;

      }

    }

  }

  // number of steps taken (?)
  k = iit + 1 - iit_vals.begin();

  if(verbose > 1){
    Rcout << "----------------------------------------" << std::endl;
    Rcout << std::endl << std::endl;
    Rcout << "sorted XB: " << std::endl << XB(iit_vals).t() << std::endl;
  }

  // initialize cut-point as the value of XB iit currently points to.
  iit_best = iit;

  // what happens if we don't have enough events or obs to split?
  // the first valid lower cut-point (at iit_vals(k)) is > the first
  // valid upper cutpoint (current value of n_risk). Put another way,
  // k (the number of steps taken from beginning of the XB vec)
  // will be > n_rows - p, where the difference on the RHS is
  // telling us where we are after taking p steps from the end
  // of the XB vec. Returning the infinite cp is a red flag.

  if(verbose > 0){
    Rcout << "j: " << j << std::endl;
    Rcout << "k: " << k << std::endl;
  }

  if (j > k){

    if(verbose > 0) {
      Rcout << "Could not find a cut-point for this XB" << std::endl;
    }

    return(R_PosInf);
  }

  if(verbose > 0){

    Rcout << "----- initializing log-rank test cutpoints -----" << std::endl;
    Rcout << "n potential cutpoints: " << k-j << std::endl;

  }

  // what happens if there are only 5 potential cut-points
  // but the value of n_split is > 5? We will just check out
  // the 5 valid cutpoints.

  // adjust k to indicate steps taken in the outer loop.
  k -= j;

  if(k > n_split){

    jit_vals = arma::linspace<arma::uvec>(0, k, n_split);

  } else {

    jit_vals = arma::linspace<arma::uvec>(0, k, k);

  }

  vec_temp.resize( jit_vals.size() );

  if(j == 0) jit_vals(jit_vals.size()-1)--;

  for(k = 0; k < vec_temp.size(); k++){
    vec_temp(k) = XB(*(iit_best - jit_vals(k)));
  }

  if(j == 0) jit_vals(jit_vals.size()-1)++;


  if(verbose > 0){

    Rcout << "cut-points chosen: ";

    Rcout << vec_temp.t();

    Rcout << "----------------------------------------" << std::endl <<
      std::endl << std::endl;

  }

  bool do_lrt = true;

  k = 0;
  j = 1;

  // begin outer loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for(jit = jit_vals.begin(); jit != jit_vals.end(); ++jit){


    if(verbose > 0){
      Rcout << "jit points to " << *jit << std::endl;
    }

    for( ; j < *jit; j++){
      group(*iit) = 1;
      --iit;
    }

    if(jit == jit_vals.begin() ||
       jit == jit_vals.end()-1){

      do_lrt = true;

    } else {

      if( vec_temp(k) == vec_temp(k+1) ||
          vec_temp(k) == vec_temp(0)   ||
          *jit <= 1){

        do_lrt = false;

      } else {

        while(XB(*iit) == XB(*(iit - 1))){

          group(*iit) = 1;
          --iit;
          ++j;

          if(verbose > 1){
            Rcout << "cutpoint dropped down one spot: ";
            Rcout << XB(*iit) << std::endl;
          }

        }

        do_lrt = true;

      }

    }

    ++k;

    if(do_lrt){

      n_risk=0;
      g_risk=0;

      observed=0;
      expected=0;

      V=0;

      break_loop = false;

      i = y_node.n_rows-1;

      if(verbose > 0){
        Rcout << "sum(group==1): " << arma::sum(group) << ";  ";
        Rcout << "sum(group==1) * w_node: " << arma::sum(group % w_node);
        Rcout << std::endl;
        Rcout << "group:" << std::endl << group(iit_vals).t() << std::endl;
      }

      // begin inner loop  - - - - - - - - - - - - -  - - - - - - - - - - - - -
      for (; ;){

        temp1 = time(i);

        n_events = 0;

        for ( ; time(i) == temp1; i--) {

          n_risk += w_node(i);
          n_events += status(i) * w_node(i);
          g_risk += group(i) * w_node(i);
          observed += status(i) * group(i) * w_node(i);

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

      if(verbose > 0){

        Rcout << "-------- log-rank test results --------" << std::endl;
        Rcout << "cutpoint: " << XB(*iit)                  << std::endl;
        Rcout << "lrt stat: " << stat_current              << std::endl;
        Rcout << "---------------------------------------" << std::endl <<
          std::endl << std::endl;

      }

      if(stat_current > stat_best){
        iit_best = iit;
        stat_best = stat_current;
      }

    }
    // end outer loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

  }

  // if the log-rank test does not detect a difference at 0.05 alpha,
  // maybe it's not a good idea to split this node.

  if(stat_best < 3.841459) return(R_PosInf);

  if(verbose > 0){
    Rcout << "Best LRT stat: " << stat_best << std::endl;
  }

  // rewind iit until it is back where it was when we got the
  // best lrt stat. While rewinding iit, also reset the group
  // values so that group is as it was when we got the best
  // lrt stat.


  while(iit <= iit_best){
    group(*iit) = 0;
    ++iit;
  }

  return(XB(*iit_best));

}

// ----------------------------------------------------------------------------
// --------------------------- ostree functions -------------------------------
// ----------------------------------------------------------------------------

void ostree_size_buffer(){

  n_slots = (nodes_max_true+1) - betas.n_cols + 10;

  if(verbose > 0){
    Rcout << "---------- buffering outputs ----------" << std::endl;
    Rcout << "ncol of betas before:  " << betas.n_cols << std::endl;
    Rcout << "number of slots added: " << n_slots      << std::endl;
  }

  betas = arma::join_horiz(
    betas,
    arma::mat(betas.n_rows, n_slots)
  );

  col_indices = arma::join_horiz(
    col_indices,
    arma::umat(col_indices.n_rows, n_slots)
  );

  children_left = arma::join_vert(
    children_left,
    arma::uvec(n_slots)
  );

  cutpoints = arma::join_vert(
    cutpoints,
    arma::vec(n_slots)
  );

  if(verbose > 0){

    Rcout << "ncol of betas after:  " << betas.n_cols  << std::endl;
    Rcout << "-------------------------------------"   << std::endl;
    Rcout << std::endl << std::endl;

  }

}

// [[Rcpp::export]]
void oobag_pred_leaf(){

  // reset values
  leaf_oobag.fill(0);

  for(i = 0; i < betas.n_cols; i++){

    if(children_left(i) != 0){

      obs_in_node = arma::find(leaf_oobag == i);

      if(obs_in_node.size() > 0){

        XB = x_oobag(obs_in_node, col_indices.col(i)) * betas.col(i);

        jit = obs_in_node.begin();

        for(j = 0; j < XB.size(); j++){

          if(XB(j) <= cutpoints(i)) {

            leaf_oobag(*jit) = children_left(i);

          } else {

            leaf_oobag(*jit) = children_left(i)+1;

          }

          jit++;

        }

        if(verbose > 0){

          arma::uvec in_left = arma::find(leaf_oobag == children_left(i));
          arma::uvec in_right = arma::find(leaf_oobag == children_left(i)+1);

          Rcout << "N to node_" << children_left(i) << ": ";
          Rcout << in_left.size() << "; ";
          Rcout << "N to node_" << children_left(i)+1 << ": ";
          Rcout << in_right.size() << std::endl;

        }

      }

    }

  }

}

// [[Rcpp::export]]
void oobag_pred_surv(){

  // allocate memory for output
  // surv_oobag.zeros(x_oobag.n_rows);

  iit_vals = arma::sort_index(leaf_oobag, "ascend");
  iit = iit_vals.begin();
  j = 0;
  k = 0;

  do {

    person_leaf = leaf_oobag(*iit);
    leaf_surv   = leaf_nodes.rows(j, leaf_node_index(k)-1);

    j = leaf_node_index(k);
    k++;

    if(verbose > 0){
      Rcout << "leaf_surv:" << std::endl << leaf_surv << std::endl;
    }

    i = 0;

    if(time_oobag < leaf_surv(leaf_surv.n_rows - 1, 0)){

      for(; i < leaf_surv.n_rows; i++){
        if (leaf_surv(i, 0) > time_oobag){
          if(i == 0)
            temp1 = 1;
          else
            temp1 = leaf_surv(i-1, 1);
          break;
        } else if (leaf_surv(i, 0) == time_oobag){
          temp1 = leaf_surv(i, 1);
          break;
        }
      }

    } else {

      // go here if prediction horizon > max time in current leaf.
      temp1 = leaf_surv(leaf_surv.n_rows - 1, 1);

    }

    // running mean: mean_k = mean_{k-1} + (new val - old val) / k
    // compute new val - old val
    // be careful, every oob row has a different denom!
    temp2 = temp1 - surv_oobag(rows_oobag(*iit));
    surv_oobag(rows_oobag(*iit)) += temp2 / denom_oobag(rows_oobag(*iit));
    ++iit;

    if(iit < iit_vals.end()){

      while(person_leaf == leaf_oobag(*iit)){

        temp2 = temp1 - surv_oobag(rows_oobag(*iit));
        surv_oobag(rows_oobag(*iit)) += temp2 / denom_oobag(rows_oobag(*iit));

        ++iit;

        if (iit == iit_vals.end()) break;

      }

    }

  } while (iit < iit_vals.end());

  if(verbose > 0){
    Rcout << "rows_oobag:" << std::endl << rows_oobag.t() << std::endl;
    Rcout << "surv_oobag:" << std::endl << surv_oobag.t() << std::endl;
  }

}


// [[Rcpp::export]]
List ostree_fit(){

  betas.fill(0);
  col_indices.fill(0);
  cutpoints.fill(0);
  children_left.fill(0);
  node_assignments.fill(0);
  leaf_nodes.fill(0);

  node_assignments.zeros(x_inbag.n_rows);
  nodes_to_grow.zeros(1);
  nodes_max_true = 0;
  leaf_node_counter = 0;
  leaf_node_index_counter = 0;


  if(verbose > 0){

    Rcout << "----------- nodes to grow -----------" << std::endl;
    Rcout << "nodes: "<< nodes_to_grow.t()           << std::endl;
    Rcout << "-------------------------------------" << std::endl <<
      std::endl << std::endl;


  }

  // ----------------------
  // ---- main do loop ----
  // ----------------------

  do {

    rows_node_combined.set_size(0);

    for(node = nodes_to_grow.begin(); node != nodes_to_grow.end(); ++node){

      if(*node >= betas.n_cols) ostree_size_buffer();

      if(nodes_to_grow[0] == 0){

        // when growing the first node, there is no need to find
        // which rows are in the node.
        rows_node = arma::linspace<arma::uvec>(0,
                                               x_inbag.n_rows-1,
                                               x_inbag.n_rows);

      } else {

        // identify which rows are in the current node.
        rows_node = arma::find(node_assignments == *node);

      }

      y_node = y_inbag.rows(rows_node);
      w_node = w_inbag(rows_node);

      if(verbose > 0){

        n_risk = arma::sum(w_node);
        n_events = arma::sum(y_node.col(1) % w_node);
        Rcout << "-------- Growing node " << *node << " --------" << std::endl;
        Rcout << "No. of observations in node: " << n_risk        << std::endl;
        Rcout << "No. of events in node:       " << n_events      << std::endl;
        Rcout << "No. of rows in node:         " << w_node.size() << std::endl;
        Rcout << "--------------------------------"               << std::endl;
        Rcout << std::endl << std::endl;

      }

      // initialize an impossible cut-point value
      // if cutpoint is still infinite later, node should not be split
      cutpoint = R_PosInf;

      // ------------------------------------------------------------------
      // ---- sample a random subset of columns with non-zero variance ----
      // ------------------------------------------------------------------

      mtry_int = mtry;
      cols_to_sample_01.fill(0);

      // constant columns are constant in the rows where events occurred

      for(j = 0; j < cols_to_sample_01.size(); j++){

        temp1 = R_PosInf;

        for(iit = rows_node.begin()+1; iit != rows_node.end(); ++iit){

          if(y_inbag.at(*iit, 1) == 1){

            if (temp1 < R_PosInf){

              if(x_inbag.at(*iit, j) != temp1){

                cols_to_sample_01[j] = 1;
                break;

              }

            } else {

              temp1 = x_inbag.at(*iit, j);

            }

          }

        }

      }

      n_cols_to_sample = arma::sum(cols_to_sample_01);

      if(n_cols_to_sample > 1){

        n_events = arma::sum(y_node.col(1));

        if(n_cols_to_sample < mtry){

          mtry_int = n_cols_to_sample;

          if(verbose > 0){
            Rcout << " ---- >=1 constant column in node rows ----" << std::endl;
            Rcout << "mtry reduced to " << mtry_temp << " from " << mtry;
            Rcout << std::endl;
            Rcout << "-------------------------------------------" << std::endl;
            Rcout << std::endl << std::endl;
          }

        }

        while(n_events / mtry_int < 2 && mtry_int > 1){
          --mtry_int;
        }

        n_cols_to_sample = mtry_int;

        if(verbose > 0){
          Rcout << "n_events: " << n_events << std::endl;
          Rcout << "mtry: " << mtry_int << std::endl;
          Rcout << "n_events per column: " << n_events/mtry_int << std::endl;
        }

        if(mtry_int > 1){

          cols_to_sample = arma::find(cols_to_sample_01);

          cols_node = Rcpp::RcppArmadillo::sample(cols_to_sample,
                                                  mtry_int,
                                                  false);

          x_node = x_inbag(rows_node, cols_node);

          n_vars = x_node.n_cols;

          if(!scale_fast){

            x_node_scale();

            if(verbose > 0){

              arma::uword temp_uword_1 = min(arma::uvec {x_node.n_rows, 5});
              Rcout << "x node scaled: " << std::endl;
              Rcout << x_node.submat(0, 0, temp_uword_1-1, x_node.n_cols-1);
              Rcout << std::endl;

            }

          }

          beta_cph = newtraph_cph();

          if(!scale_fast){
            for(i = 0; i < x_transforms.n_rows; i++){
              x_node.col(i) /= x_transforms(i,1);
              x_node.col(i) += x_transforms(i,0);
            }
          }

          if(arma::any(beta_cph)){

            if(verbose > 0){

              arma::uword temp_uword_1 = min(arma::uvec {x_node.n_rows, 5});
              Rcout << "x node unscaled: " << std::endl;
              Rcout << x_node.submat(0, 0, temp_uword_1-1, x_node.n_cols-1);
              Rcout << std::endl;

            }

            XB = x_node * beta_cph;
            cutpoint = lrt_multi();

          }

        }

      }

      if(!std::isinf(cutpoint)){

        rows_node_combined = arma::join_cols(rows_node_combined,
                                             rows_node);

        nn_left   = nodes_max_true + 1;
        nodes_max_true = nodes_max_true + 2;

        if(verbose > 0){

          Rcout << "-------- New nodes created --------" << std::endl;
          Rcout << "Left node: " << nn_left              << std::endl;
          Rcout << "Right node: " << nodes_max_true      << std::endl;
          Rcout << "-----------------------------------" << std::endl <<
            std::endl << std::endl;

        }

        i = 0;

        for(iit = rows_node.begin(); iit != rows_node.end(); ++iit){

          node_assignments[*iit] = nn_left + group[i];

          i++;

        }

        for(i = 0; i < n_cols_to_sample; i++){
          betas.at(i, *node) = beta_cph[i];
          col_indices.at(i, *node) = cols_node(i);
        }

        children_left[*node] = nn_left;
        cutpoints[*node] = cutpoint;

      } else {

        leaf_surv_small(y_node, w_node);

        if(verbose > 0){
          Rcout << "-------- creating a new leaf --------" << std::endl;
          Rcout << "n_obs:    "  << arma::sum(w_node)      << std::endl;
          Rcout << "n_events: "  << arma::sum(w_node % y_node.col(1));
          Rcout << std::endl;
          Rcout << "------------------------------------"  << std::endl;
          Rcout << std::endl << std::endl;
        }

        // leaf_nodes[make_node_name(*node)] = leaf;

      }

    }

    node_sums.zeros(nodes_max_true, 2);
    iit = rows_node_combined.begin();


    for(; iit < rows_node_combined.end(); ++iit){
      node_sums(node_assignments(*iit)-1, 0) += y_inbag(*iit,1) * w_inbag(*iit);
      node_sums(node_assignments(*iit)-1, 1) += w_inbag(*iit);
    }

    if(verbose > 0){

      vec_temp = arma::regspace<arma::vec>(1, 1, node_sums.n_rows);

      arma::umat node_sums_temp = arma::conv_to<arma::umat>::from(
        arma::join_horiz(vec_temp, node_sums)
      );

      Rcout << "node_sums: " << std::endl;
      Rcout << arma::round(node_sums_temp) << std::endl;

    }

    nodes_to_grow.zeros(nodes_max_true);

    for(i = 0; i < node_sums.n_rows; i++){

      if(node_sums(i, 0) >= 2 * leaf_min_events &&
         node_sums(i, 1) >= 2 * leaf_min_obs){
        // split it

        // use i+1; the ith row of node_sums is the i+1 node
        nodes_to_grow(i) = i + 1;

      } else if (node_sums(i, 0) > 0 && node_sums(i, 1) > 0) {
        // a new leaf

        // use i+1; nodes starts at 1 and i starts at 0
        rows_leaf    = arma::find(node_assignments == i+1);

        leaf_surv_small(y_inbag.rows(rows_leaf),
                        w_inbag(rows_leaf));


        if(verbose > 0){
          Rcout << "-------- creating a new leaf --------" << std::endl;
          Rcout << "name: node_" << i+1                    << std::endl;
          Rcout << "n_obs:    "  << arma::sum(w_inbag(rows_leaf));
          Rcout << std::endl;
          Rcout << "n_events: ";
          vec_temp = y_inbag.col(1);
          Rcout << arma::sum(w_inbag(rows_leaf) % vec_temp(rows_leaf));
          Rcout << std::endl;
          Rcout << "------------------------------------";
          Rcout << std::endl << std::endl << std::endl;
        }

        // leaf_nodes[make_node_name(i+1)] = leaf;

      }

    }

    nodes_to_grow = nodes_to_grow(arma::find(nodes_to_grow));

  } while (nodes_to_grow.size() > 0);

  return(
    List::create(
      _["leaf_nodes"] = leaf_nodes.rows(arma::span(0, leaf_node_counter-1)),
      _["leaf_node_index"] = leaf_node_index.rows(arma::span(0, leaf_node_index_counter-1)),
      _["betas"] = betas.cols(arma::span(0, nodes_max_true)),
      _["col_indices"] = col_indices.cols(arma::span(0, nodes_max_true)),
      _["cut_points"] = cutpoints(arma::span(0, nodes_max_true)),
      _["children_left"] = children_left(arma::span(0, nodes_max_true))
    )
  );


}

// [[Rcpp::export]]
List orsf_fit(NumericMatrix&   x,
              NumericMatrix&   y,
              const int&       n_tree = 2,
              const int&       n_split_ = 5,
              const int&       mtry_ = 4,
              const double&    leaf_min_events_ = 5,
              const double&    leaf_min_obs_ = 10,
              const int&       cph_method_ = 1,
              const double&    cph_eps_ = 1e-8,
              const int&       cph_iter_max_ = 7,
              const double&    cph_pval_max_ = 0.95,
              const bool&      oobag_pred_ = false,
              const bool&      scale_fast_ = true){


  // convert inputs into arma objects
  x_input = arma::mat(x.begin(), x.nrow(), x.ncol(), false);
  y_input = arma::mat(y.begin(), y.nrow(), y.ncol(), false);

  // these change later in ostree_fit()
  n_rows = x_input.n_rows;
  n_vars = x_input.n_cols;

  if(verbose > 0){
    Rcout << "------------ dimensions ------------"  << std::endl;
    Rcout << "N obs total: "     << n_rows           << std::endl;
    Rcout << "N columns total: " << n_vars           << std::endl;
    Rcout << "------------------------------------";
    Rcout << std::endl << std::endl << std::endl;
  }

  n_split          = n_split_;
  mtry             = mtry_;
  leaf_min_events  = leaf_min_events_;
  leaf_min_obs     = leaf_min_obs_;
  cph_method       = cph_method_;
  cph_eps          = cph_eps_;
  cph_iter_max     = cph_iter_max_;
  cph_pval_max     = cph_pval_max_;
  oobag_pred       = oobag_pred_;
  scale_fast       = scale_fast_;
  temp1            = 1.0 / n_rows;

  if(oobag_pred){ time_oobag = arma::median(y_input.col(0)); }

  if(verbose > 0){
    Rcout << "------------ input variables ------------" << std::endl;
    Rcout << "n_split: "         << n_split              << std::endl;
    Rcout << "mtry: "            << mtry                 << std::endl;
    Rcout << "leaf_min_events: " << leaf_min_events      << std::endl;
    Rcout << "leaf_min_obs: "    << leaf_min_obs         << std::endl;
    Rcout << "cph_method: "      << cph_method           << std::endl;
    Rcout << "cph_eps: "         << cph_eps              << std::endl;
    Rcout << "cph_iter_max: "    << cph_iter_max         << std::endl;
    Rcout << "cph_pval_max: "    << cph_pval_max         << std::endl;
    Rcout << "-----------------------------------------" << std::endl;
    Rcout << std::endl << std::endl;
  }

  if(scale_fast){ x_input_scale(); }

  // ----------------------------------------------------
  // ---- sample weights to mimic a bootstrap sample ----
  // ----------------------------------------------------

  // s is the number of times you might get selected into
  // a bootstrap sample. Realistically this won't be >10,
  // but it could technically be as big as n_row.
  IntegerVector s = seq(0, 10);

  // compute probability of being selected into the bootstrap
  // 0 times, 1, times, ..., 9 times, or 10 times.
  NumericVector probs = dbinom(s, n_rows, temp1, false);

  // ---------------------------------------------
  // ---- preallocate memory for tree outputs ----
  // ---------------------------------------------



  cols_to_sample_01.zeros(n_vars);
  leaf_nodes.zeros(n_rows, 2);

  if(oobag_pred){
    surv_oobag.zeros(n_rows);
    denom_oobag.zeros(n_rows);
  }

  // guessing the number of nodes needed to grow a tree
  nodes_max_guess = std::ceil(n_rows / leaf_min_events);

  betas.zeros(mtry, nodes_max_guess);
  col_indices.zeros(mtry, nodes_max_guess);
  cutpoints.zeros(nodes_max_guess);
  children_left.zeros(nodes_max_guess);
  leaf_node_index.zeros(nodes_max_guess);


  List forest(n_tree);

  for(int tree = 0; tree < n_tree; tree++){

    // --------------------------------------------
    // ---- initialize parameters to grow tree ----
    // --------------------------------------------

    w_inbag    = as<arma::vec>(sample(s, n_rows, true, probs));
    rows_inbag = arma::find(w_inbag != 0);
    rows_oobag = arma::find(w_inbag == 0);
    w_inbag    = w_inbag(rows_inbag);

    if(verbose > 0){

      Rcout << "------------ boot weights ------------" << std::endl;
      Rcout << "pr(inbag): " << 1-pow(1-temp1,n_rows)   << std::endl;
      Rcout << "total: "     << arma::sum(w_inbag)      << std::endl;
      Rcout << "N > 0: "     << rows_inbag.size()       << std::endl;
      Rcout << "--------------------------------------" <<
        std::endl << std::endl << std::endl;

    }

    x_inbag = x_input.rows(rows_inbag);
    y_inbag = y_input.rows(rows_inbag);

    if(oobag_pred){
      x_oobag = x_input.rows(rows_oobag);
      y_oobag = y_input.rows(rows_oobag);
      leaf_oobag.set_size(rows_oobag.size());
    }

    if(verbose > 0){

      arma::uword temp_uword_1, temp_uword_2;

      if(x_inbag.n_rows < 5)
        temp_uword_1 = x_inbag.n_rows-1;
      else
        temp_uword_1 = 5;

      if(x_inbag.n_cols < 5)
        temp_uword_2 = x_inbag.n_cols-1;
      else
        temp_uword_2 = 5;

      Rcout << "x inbag: " << std::endl <<
        x_inbag.submat(0, 0,
                       temp_uword_1,
                       temp_uword_2) << std::endl;

    }

    forest[tree] = ostree_fit();

    if(oobag_pred){

      denom_oobag(rows_oobag) += 1;
      oobag_pred_leaf();
      oobag_pred_surv();

    }

  }

  return(
    List::create(
      _["forest"] = forest,
      _["surv_oobag"] = surv_oobag,
      _["scale_fast"] = scale_fast,
      _["mtry"] = mtry
    )
  );


}



