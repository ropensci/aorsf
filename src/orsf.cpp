#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;

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
 // the n_ variables could be integers but it
 // is safer and faster when they are doubles
 n_events,
 n_events_total,
 n_events_right,
 n_events_left,
 n_risk,
 n_risk_right,
 n_risk_left,
 n_risk_sub,
 g_risk,
 temp1,
 temp2,
 halving,
 stat_current,
 stat_best,
 w_node_person,
 xb,
 risk,
 loglik,
 cutpoint,
 observed,
 expected,
 V,
 leaf_min_obs,
 leaf_min_events,
 split_min_events,
 split_min_obs,
 split_min_stat,
 time_pred,
 ll_second,
 ll_init,
 net_alpha;

int
 // verbose=0,
 max_retry,
 n_retry,
 tree,
 mtry_int,
 net_df_target,
 oobag_eval_every;

char type_beta, type_oobag_eval, oobag_importance_type;

// armadillo unsigned integers
uword
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
 cph_method,
 cph_iter_max,
 n_split,
 nodes_max_guess,
 nodes_max_true,
 n_cols_to_sample,
 nn_left,
 leaf_node_counter,
 leaf_node_index_counter,
 oobag_eval_counter;

bool
 break_loop, // a delayed break statement
 oobag_pred,
 oobag_importance,
 use_tree_seed,
 cph_do_scale;

// armadillo vectors (doubles)
vec
 vec_temp,
 times_pred,
 eval_oobag,
 node_assignments,
 nodes_grown,
 surv_pvec,
 surv_pvec_output,
 denom_pred,
 beta_current,
 beta_new,
 beta_fit,
 vi_pval_numer,
 vi_pval_denom,
 cutpoints,
 w_input,
 w_inbag,
 w_node,
 group,
 u,
 a,
 a2,
 XB,
 Risk;

// armadillo unsigned integer vectors
uvec
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
 nodes_to_grow_next,
 obs_in_node,
 children_left,
 leaf_pred;

// armadillo iterators for unsigned integer vectors
uvec::iterator
 iit,
 iit_best,
 jit,
 node;

// armadillo matrices (doubles)
mat
 x_input,
 x_transforms,
 y_input,
 x_inbag,
 y_inbag,
 x_node,
 y_node,
 x_pred,
 // x_mean,
 vmat,
 cmat,
 cmat2,
 betas,
 leaf_node,
 leaf_nodes,
 surv_pmat;

umat
 col_indices,
 leaf_indices;

cube
 surv_pcube;

List ostree;

NumericMatrix
 beta_placeholder,
 xx,
 yy;

CharacterVector yy_names = CharacterVector::create("time","status");

NumericVector ww;

Environment base_env("package:base");

Function set_seed_r = base_env["set.seed"];

// [[Rcpp::export]]
arma::uvec std_setdiff(arma::uvec& x, arma::uvec& y) {

 std::vector<int> a = conv_to< std::vector<int> >::from(sort(x));
 std::vector<int> b = conv_to< std::vector<int> >::from(sort(y));
 std::vector<int> out;

 std::set_difference(a.begin(), a.end(),
                     b.begin(), b.end(),
                     std::inserter(out, out.end()));

 return conv_to<uvec>::from(out);

}

// ----------------------------------------------------------------------------
// ---------------------------- scaling functions -----------------------------
// ----------------------------------------------------------------------------

void x_node_scale(){

 // set aside memory for outputs
 // first column holds the mean values
 // second column holds the scale values

 x_transforms.zeros(n_vars, 2);
 vec means  = x_transforms.unsafe_col(0);   // Reference to column 1
 vec scales = x_transforms.unsafe_col(1);   // Reference to column 2

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

void x_node_means(){

 x_transforms.zeros(n_vars, 1);
 w_node_sum = sum(w_node);

 for(i = 0; i < n_vars; i++) {

  x_transforms.at(i, 0) = sum( w_node % x_node.col(i) ) / w_node_sum;

 }

}

// [[Rcpp::export]]
List x_node_scale_exported(NumericMatrix& x_,
                           NumericVector& w_){

 x_node = mat(x_.begin(), x_.nrow(), x_.ncol(), false);
 w_node = vec(w_.begin(), w_.length(), false);
 n_vars = x_node.n_cols;

 x_node_scale();

 return(
  List::create(
   _["x_scaled"] = x_node,
   _["x_transforms"] = x_transforms
  )
 );

}

// ----------------------------------------------------------------------------
// -------------------------- leaf_surv functions -----------------------------
// ----------------------------------------------------------------------------

void leaf_kaplan(const arma::mat& y,
                 const arma::vec& w){

 leaf_indices(leaf_node_index_counter, 1) = leaf_node_counter;
 i = leaf_node_counter;

 // find the first unique event time
 person = 0;

 while(y.at(person, 1) == 0){
  person++;
 }

 // now person corresponds to the first event time
 leaf_nodes.at(i, 0) = y.at(person, 0);  // see above
 temp2 = y.at(person, 0);

 i++;

 // find the rest of the unique event times
 for( ; person < y.n_rows; person++){

  if(temp2 != y.at(person, 0) && y.at(person, 1) == 1){

   leaf_nodes.at(i, 0) = y.at(person,0);
   temp2 = y.at(person, 0);
   i++;

  }

 }


 // reset for kaplan meier loop
 n_risk = sum(w);
 person = 0;
 temp1 = 1.0;

 do {

  n_events   = 0;
  n_risk_sub = 0;
  temp2      = y.at(person, 0);

  while(y.at(person, 0) == temp2){

   n_risk_sub += w.at(person);
   n_events += y.at(person, 1) * w.at(person);

   if(person == y.n_rows-1) break;

   person++;

  }

  // only do km if a death was observed

  if(n_events > 0){

   temp1 = temp1 * (n_risk - n_events) / n_risk;
   leaf_nodes.at(leaf_node_counter, 1) = temp1;
   leaf_node_counter++;

  }

  n_risk -= n_risk_sub;

 } while (leaf_node_counter < i);


 leaf_indices(leaf_node_index_counter, 2) = leaf_node_counter-1;
 leaf_node_index_counter++;

 if(leaf_node_index_counter >= leaf_indices.n_rows){
  leaf_indices.insert_rows(leaf_indices.n_rows, 10);
 }

}

// [[Rcpp::export]]
arma::mat leaf_kaplan_testthat(const arma::mat& y,
                               const arma::vec& w){


 leaf_nodes.set_size(y.n_rows, 2);
 leaf_node_counter = 0;

 // find the first unique event time
 person = 0;

 while(y.at(person, 1) == 0){
  person++;
 }

 // now person corresponds to the first event time
 leaf_nodes.at(leaf_node_counter, 0) = y.at(person, 0);  // see above
 temp2 = y.at(person, 0);

 leaf_node_counter++;

 // find the rest of the unique event times
 for( ; person < y.n_rows; person++){

  if(temp2 != y.at(person, 0) && y.at(person, 1) == 1){

   leaf_nodes.at(leaf_node_counter, 0) = y.at(person,0);
   temp2 = y.at(person, 0);
   leaf_node_counter++;

  }

 }


 // reset for kaplan meier loop
 i = leaf_node_counter;
 n_risk = sum(w);
 person = 0;
 temp1 = 1.0;
 leaf_node_counter = 0;


 do {

  n_events   = 0;
  n_risk_sub = 0;
  temp2      = y.at(person, 0);

  while(y.at(person, 0) == temp2){

   n_risk_sub += w.at(person);
   n_events += y.at(person, 1) * w.at(person);

   if(person == y.n_rows-1) break;

   person++;

  }

  // only do km if a death was observed

  if(n_events > 0){

   temp1 = temp1 * (n_risk - n_events) / n_risk;
   leaf_nodes.at(leaf_node_counter, 1) = temp1;
   leaf_node_counter++;

  }

  n_risk -= n_risk_sub;

 } while (leaf_node_counter < i);

 leaf_nodes.resize(leaf_node_counter, 2);

 return(leaf_nodes);

}




// ----------------------------------------------------------------------------
// ---------------------------- cholesky functions ----------------------------
// ----------------------------------------------------------------------------

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

  if (pivot < R_PosInf && pivot > eps_chol) {

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

void cholesky_invert(){

 /*
  ** invert the cholesky in the lower triangle
  **   take full advantage of the cholesky's diagonal of 1's
  */
 for (i=0; i<n_vars; i++){

  if (vmat.at(i,i) >0) {

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

double newtraph_cph_init(){

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

arma::vec newtraph_cph(){

 beta_current.zeros(n_vars);
 beta_new.zeros(n_vars);

 // these are filled with initial values later
 XB.set_size(x_node.n_rows);
 Risk.set_size(x_node.n_rows);
 u.set_size(n_vars);
 a.set_size(n_vars);
 a2.set_size(n_vars);
 vmat.set_size(n_vars, n_vars);
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

  if(oobag_importance_type == 'A'){

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


// [[Rcpp::export]]
arma::vec newtraph_cph_testthat(NumericMatrix& x_in,
                                NumericMatrix& y_in,
                                NumericVector& w_in,
                                int method,
                                double cph_eps_,
                                int iter_max){


 x_node = mat(x_in.begin(), x_in.nrow(), x_in.ncol(), false);
 y_node = mat(y_in.begin(), y_in.nrow(), y_in.ncol(), false);
 w_node = vec(w_in.begin(), w_in.length(), false);

 cph_do_scale = true;

 cph_method = method;
 cph_eps = cph_eps_;
 cph_iter_max = iter_max;
 n_vars = x_node.n_cols;

 vi_pval_numer.zeros(x_node.n_cols);
 vi_pval_denom.zeros(x_node.n_cols);
 cols_node = regspace<uvec>(0, x_node.n_cols - 1);

 x_node_scale();

 vec out = newtraph_cph();

 return(out);

}

// ----------------------------------------------------------------------------
// ---------------------------- node functions --------------------------------
// ----------------------------------------------------------------------------

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
 iit_vals = sort_index(XB, "ascend");

 // unsafe columns point to cols in y_node.
 vec y_status = y_node.unsafe_col(1);
 vec y_time = y_node.unsafe_col(0);

 // first determine the lowest value of XB that will
 // be a valid cut-point to split a node. A valid cut-point
 // is one that, if used, will result in at least leaf_min_obs
 // and leaf_min_events in both the left and right node.

 n_events = 0;
 n_risk = 0;

 // if(verbose > 1){
 //  Rcout << "----- finding cut-point boundaries -----" << std::endl;
 // }

 // Iterate through the sorted values of XB, in ascending order.

 for(iit = iit_vals.begin(); iit < iit_vals.end()-1; ++iit){

  n_events += y_status[*iit] * w_node[*iit];
  n_risk += w_node[*iit];

  // If we want to make the current value of XB a cut-point, we need
  // to make sure the next value of XB isn't equal to this current value.
  // Otherwise, we will have the same value of XB in both groups!

  // if(verbose > 1){
  //  Rcout << XB[*iit]     << " ---- ";
  //  Rcout << XB[*(iit+1)] << " ---- ";
  //  Rcout << n_events     << " ---- ";
  //  Rcout << n_risk       << std::endl;
  // }

  if(XB[*iit] != XB[*(iit+1)]){

   // if(verbose > 1){
   //  Rcout << "********* New cut-point here ********" << std::endl;
   // }


   if( n_events >= leaf_min_events &&
       n_risk   >= leaf_min_obs) {

    // if(verbose > 1){
    //  Rcout << std::endl;
    //  Rcout << "lower cutpoint: "         << XB[*iit] << std::endl;
    //  Rcout << " - n_events, left node: " << n_events << std::endl;
    //  Rcout << " - n_risk, left node:   " << n_risk   << std::endl;
    //  Rcout << std::endl;
    // }

    break;

   }

  }

 }

 // if(verbose > 1){
 //  if(iit >= iit_vals.end()-1) {
 //   Rcout << "Could not find a valid lower cut-point" << std::endl;
 //  }
 // }


 j = iit - iit_vals.begin();

 // got to reset these before finding the upper limit
 n_events=0;
 n_risk=0;

 // do the first step in the loop manually since we need to
 // refer to iit+1 in all proceeding steps.

 for(iit = iit_vals.end()-1; iit >= iit_vals.begin()+1; --iit){

  n_events += y_status[*iit] * w_node[*iit];
  n_risk   += w_node[*iit];
  group[*iit] = 1;

  // if(verbose > 1){
  //  Rcout << XB[*iit]     << " ---- ";
  //  Rcout << XB(*(iit-1)) << " ---- ";
  //  Rcout << n_events     << " ---- ";
  //  Rcout << n_risk       << std::endl;
  // }

  if ( XB[*iit] != XB[*(iit-1)] ) {

   // if(verbose > 1){
   //  Rcout << "********* New cut-point here ********" << std::endl;
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

    // if(verbose > 1){
    //  Rcout << std::endl;
    //  Rcout << "upper cutpoint: " << XB[*iit] << std::endl;
    //  Rcout << " - n_events, right node: " << n_events    << std::endl;
    //  Rcout << " - n_risk, right node:   " << n_risk      << std::endl;
    // }

    break;

   }

  }

 }

 // number of steps taken
 k = iit + 1 - iit_vals.begin();

 // if(verbose > 1){
 //  Rcout << "----------------------------------------" << std::endl;
 //  Rcout << std::endl << std::endl;
 //  Rcout << "sorted XB: " << std::endl << XB(iit_vals).t() << std::endl;
 // }

 // initialize cut-point as the value of XB iit currently points to.
 iit_best = iit;

 // what happens if we don't have enough events or obs to split?
 // the first valid lower cut-point (at iit_vals(k)) is > the first
 // valid upper cutpoint (current value of n_risk). Put another way,
 // k (the number of steps taken from beginning of the XB vec)
 // will be > n_rows - p, where the difference on the RHS is
 // telling us where we are after taking p steps from the end
 // of the XB vec. Returning the infinite cp is a red flag.

 // if(verbose > 1){
 //  Rcout << "j: " << j << std::endl;
 //  Rcout << "k: " << k << std::endl;
 // }

 if (j > k){

  // if(verbose > 1) {
  //  Rcout << "Could not find a cut-point for this XB" << std::endl;
  // }

  return(R_PosInf);
 }

 // if(verbose > 1){
 //
 //  Rcout << "----- initializing log-rank test cutpoints -----" << std::endl;
 //  Rcout << "n potential cutpoints: " << k-j << std::endl;
 //
 // }


 // adjust k to indicate the number of valid cut-points
 k -= j;

 if(k > n_split){

  jit_vals = linspace<uvec>(0, k, n_split);

 } else {

  // what happens if there are only 5 potential cut-points
  // but the value of n_split is > 5? We will just check out
  // the 5 valid cutpoints.
  jit_vals = linspace<uvec>(0, k, k);

 }

 vec_temp.resize( jit_vals.size() );

 // protection from going out of bounds with jit_vals(k) below
 if(j == 0) jit_vals(jit_vals.size()-1)--;

 // put the indices of potential cut-points into vec_temp
 for(k = 0; k < vec_temp.size(); k++){
  vec_temp[k] = XB(*(iit_best - jit_vals[k]));
 }

 // back to how it was!
 if(j == 0) jit_vals(jit_vals.size()-1)++;

 // if(verbose > 1){
 //
 //  Rcout << "cut-points chosen: ";
 //
 //  Rcout << vec_temp.t();
 //
 //  Rcout << "----------------------------------------" << std::endl <<
 //   std::endl << std::endl;
 //
 // }

 bool do_lrt = true;

 k = 0;
 j = 1;

 // begin outer loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for(jit = jit_vals.begin(); jit != jit_vals.end(); ++jit){


  // if(verbose > 1){
  //  Rcout << "jit points to " << *jit << std::endl;
  // }

  // switch group values from 0 to 1 until you get to the next cut-point
  for( ; j < *jit; j++){
   group[*iit] = 1;
   --iit;
  }

  if(jit == jit_vals.begin() ||
     jit == jit_vals.end()-1){

   do_lrt = true;

  } else {

   if( vec_temp[k] == vec_temp[k+1] ||
       vec_temp[k] == vec_temp[0]   ||
       *jit <= 1){

       do_lrt = false;

   } else {

    while( XB[*iit] == XB[*(iit - 1)] ){

     group[*iit] = 1;
     --iit;
     ++j;

     // if(verbose > 1){
     //  Rcout << "cutpoint dropped down one spot: ";
     //  Rcout << XB[*iit] << std::endl;
     // }

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

   // if(verbose > 1){
   //  Rcout << "sum(group==1): " << sum(group) << ";  ";
   //  Rcout << "sum(group==1 * w_node): " << sum(group % w_node);
   //  Rcout << std::endl;
   //  if(verbose > 1){
   //   Rcout << "group:" << std::endl;
   //   Rcout << group(iit_vals).t() << std::endl;
   //  }
   // }


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

   // if(verbose > 1){
   //
   //  Rcout << "-------- log-rank test results --------" << std::endl;
   //  Rcout << "cutpoint: " << XB[*iit]                  << std::endl;
   //  Rcout << "lrt stat: " << stat_current              << std::endl;
   //  Rcout << "---------------------------------------" << std::endl <<
   //   std::endl << std::endl;
   //
   // }

   if(stat_current > stat_best){
    iit_best = iit;
    stat_best = stat_current;
    n_events_right = observed;
    n_risk_right = g_risk;
    n_risk_left = n_risk - g_risk;
   }

  }
  // end outer loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

 }

 // if the log-rank test does not detect a difference at 0.05 alpha,
 // maybe it's not a good idea to split this node.

 if(stat_best < split_min_stat) return(R_PosInf);

 // if(verbose > 1){
 //  Rcout << "Best LRT stat: " << stat_best << std::endl;
 // }

 // rewind iit until it is back where it was when we got the
 // best lrt stat. While rewinding iit, also reset the group
 // values so that group is as it was when we got the best
 // lrt stat.


 while(iit <= iit_best){
  group[*iit] = 0;
  ++iit;
 }

 return(XB[*iit_best]);

}

// [[Rcpp::export]]
List lrt_multi_testthat(NumericMatrix& y_node_,
                        NumericVector& w_node_,
                        NumericVector& XB_,
                        int n_split_,
                        int leaf_min_events_,
                        int leaf_min_obs_
){

 y_node = mat(y_node_.begin(), y_node_.nrow(), y_node_.ncol(), false);
 w_node = vec(w_node_.begin(), w_node_.length(), false);
 XB = vec(XB_.begin(), XB_.length(), false);

 n_split = n_split_;
 leaf_min_events = leaf_min_events_;
 leaf_min_obs = leaf_min_obs_;

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

 vec cutpoints_used(n_split);
 vec lrt_statistics(n_split);
 uword list_counter = 0;

 // group should be initialized as all 0s
 group.zeros(y_node.n_rows);

 // initialize at the lowest possible LRT stat value
 stat_best = 0;

 // sort XB- we need to iterate over the sorted indices
 iit_vals = sort_index(XB, "ascend");

 // unsafe columns point to cols in y_node.
 vec y_status = y_node.unsafe_col(1);
 vec y_time = y_node.unsafe_col(0);

 // first determine the lowest value of XB that will
 // be a valid cut-point to split a node. A valid cut-point
 // is one that, if used, will result in at least leaf_min_obs
 // and leaf_min_events in both the left and right node.

 n_events = 0;
 n_risk = 0;

 // if(verbose > 1){
 //  Rcout << "----- finding cut-point boundaries -----" << std::endl;
 // }

 // Iterate through the sorted values of XB, in ascending order.

 for(iit = iit_vals.begin(); iit < iit_vals.end()-1; ++iit){

  n_events += y_status(*iit) * w_node(*iit);
  n_risk += w_node(*iit);

  // If we want to make the current value of XB a cut-point, we need
  // to make sure the next value of XB isn't equal to this current value.
  // Otherwise, we will have the same value of XB in both groups!

  // if(verbose > 1){
  //  Rcout << XB(*iit)     << " ---- ";
  //  Rcout << XB(*(iit+1)) << " ---- ";
  //  Rcout << n_events     << " ---- ";
  //  Rcout << n_risk       << std::endl;
  // }

  if(XB(*iit) != XB(*(iit+1))){

   // if(verbose > 1){
   //  Rcout << "********* New cut-point here ********" << std::endl;
   // }


   if( n_events >= leaf_min_events &&
       n_risk   >= leaf_min_obs) {

    // if(verbose > 1){
    //  Rcout << std::endl;
    //  Rcout << "lower cutpoint: "         << XB(*iit) << std::endl;
    //  Rcout << " - n_events, left node: " << n_events << std::endl;
    //  Rcout << " - n_risk, left node:   " << n_risk   << std::endl;
    //  Rcout << std::endl;
    // }

    break;

   }

  }

 }

 // if(verbose > 1){
 //  if(iit >= iit_vals.end()-1) {
 //   Rcout << "Could not find a valid lower cut-point" << std::endl;
 //  }
 // }


 j = iit - iit_vals.begin();

 // got to reset these before finding the upper limit
 n_events=0;
 n_risk=0;

 // do the first step in the loop manually since we need to
 // refer to iit+1 in all proceeding steps.

 for(iit = iit_vals.end()-1; iit >= iit_vals.begin()+1; --iit){

  n_events += y_status(*iit) * w_node(*iit);
  n_risk   += w_node(*iit);
  group(*iit) = 1;

  // if(verbose > 1){
  //  Rcout << XB(*iit)     << " ---- ";
  //  Rcout << XB(*(iit-1)) << " ---- ";
  //  Rcout << n_events     << " ---- ";
  //  Rcout << n_risk       << std::endl;
  // }

  if(XB(*iit) != XB(*(iit-1))){

   // if(verbose > 1){
   //  Rcout << "********* New cut-point here ********" << std::endl;
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

    // if(verbose > 1){
    //  Rcout << std::endl;
    //  Rcout << "upper cutpoint: " << XB(*iit) << std::endl;
    //  Rcout << " - n_events, right node: " << n_events    << std::endl;
    //  Rcout << " - n_risk, right node:   " << n_risk      << std::endl;
    // }

    break;

   }

  }

 }

 // number of steps taken
 k = iit + 1 - iit_vals.begin();

 // if(verbose > 1){
 //  Rcout << "----------------------------------------" << std::endl;
 //  Rcout << std::endl << std::endl;
 //  Rcout << "sorted XB: " << std::endl << XB(iit_vals).t() << std::endl;
 // }

 // initialize cut-point as the value of XB iit currently points to.
 iit_best = iit;

 // what happens if we don't have enough events or obs to split?
 // the first valid lower cut-point (at iit_vals(k)) is > the first
 // valid upper cutpoint (current value of n_risk). Put another way,
 // k (the number of steps taken from beginning of the XB vec)
 // will be > n_rows - p, where the difference on the RHS is
 // telling us where we are after taking p steps from the end
 // of the XB vec. Returning the infinite cp is a red flag.

 // if(verbose > 1){
 //  Rcout << "j: " << j << std::endl;
 //  Rcout << "k: " << k << std::endl;
 // }

 if (j > k){

  // if(verbose > 1) {
  //  Rcout << "Could not find a cut-point for this XB" << std::endl;
  // }

  return(R_PosInf);
 }

 // if(verbose > 1){
 //
 //  Rcout << "----- initializing log-rank test cutpoints -----" << std::endl;
 //  Rcout << "n potential cutpoints: " << k-j << std::endl;
 //
 // }

 // what happens if there are only 5 potential cut-points
 // but the value of n_split is > 5? We will just check out
 // the 5 valid cutpoints.

 // adjust k to indicate steps taken in the outer loop.
 k -= j;

 if(k > n_split){

  jit_vals = linspace<uvec>(0, k, n_split);

 } else {

  jit_vals = linspace<uvec>(0, k, k);

 }

 vec_temp.resize( jit_vals.size() );

 if(j == 0) jit_vals(jit_vals.size()-1)--;

 for(k = 0; k < vec_temp.size(); k++){
  vec_temp(k) = XB(*(iit_best - jit_vals(k)));
 }

 if(j == 0) jit_vals(jit_vals.size()-1)++;


 // if(verbose > 1){
 //
 //  Rcout << "cut-points chosen: ";
 //
 //  Rcout << vec_temp.t();
 //
 //  Rcout << "----------------------------------------" << std::endl <<
 //   std::endl << std::endl;
 //
 // }

 bool do_lrt = true;

 k = 0;
 j = 1;

 // begin outer loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for(jit = jit_vals.begin(); jit != jit_vals.end(); ++jit){


  // if(verbose > 1){
  //  Rcout << "jit points to " << *jit << std::endl;
  // }

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

     // if(verbose > 1){
     //  Rcout << "cutpoint dropped down one spot: ";
     //  Rcout << XB(*iit) << std::endl;
     // }

    }

    do_lrt = true;

   }

  }

  ++k;

  if(do_lrt){

   cutpoints_used(list_counter) = XB(*iit);

   n_risk=0;
   g_risk=0;

   observed=0;
   expected=0;

   V=0;

   break_loop = false;

   i = y_node.n_rows-1;

   // if(verbose > 1){
   //  Rcout << "sum(group==1): " << sum(group) << ";  ";
   //  Rcout << "sum(group==1 * w_node): " << sum(group % w_node);
   //  Rcout << std::endl;
   //  if(verbose > 1){
   //   Rcout << "group:" << std::endl;
   //   Rcout << group(iit_vals).t() << std::endl;
   //  }
   // }


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

   // if(verbose > 1){
   //
   //  Rcout << "-------- log-rank test results --------" << std::endl;
   //  Rcout << "cutpoint: " << XB(*iit)                  << std::endl;
   //  Rcout << "lrt stat: " << stat_current              << std::endl;
   //  Rcout << "---------------------------------------" << std::endl <<
   //   std::endl << std::endl;
   //
   // }

   if(stat_current > stat_best){
    iit_best = iit;
    stat_best = stat_current;
    n_events_right = observed;
    n_risk_right = g_risk;
    n_risk_left = n_risk - g_risk;
   }

  }
  // end outer loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

 }

 // if the log-rank test does not detect a difference at 0.05 alpha,
 // maybe it's not a good idea to split this node.

 if(stat_best < 3.841459) return(R_PosInf);

 // if(verbose > 1){
 //  Rcout << "Best LRT stat: " << stat_best << std::endl;
 // }

 // rewind iit until it is back where it was when we got the
 // best lrt stat. While rewinding iit, also reset the group
 // values so that group is as it was when we got the best
 // lrt stat.


 while(iit <= iit_best){
  group(*iit) = 0;
  ++iit;
 }

 return(List::create(_["cutpoints"] = cutpoints_used,
                     _["statistic"] = lrt_statistics));

}



void oobag_pred_surv_uni(){

 iit_vals = sort_index(leaf_pred, "ascend");
 iit = iit_vals.begin();

 do {

  person_leaf = leaf_pred[*iit];

  // find the current leaf
  for(i = 0; i < leaf_indices.n_rows; i++){
   if(leaf_indices.at(i, 0) == person_leaf){
    break;
   }
  }

  // get submat view for this leaf
  leaf_node = leaf_nodes.rows(leaf_indices(i, 1),
                              leaf_indices(i, 2));

  // if(verbose > 1){
  //  Rcout << "leaf_node:" << std::endl << leaf_node << std::endl;
  // }

  i = 0;

  if(time_pred < leaf_node.at(leaf_node.n_rows - 1, 0)){

   for(; i < leaf_node.n_rows; i++){
    if (leaf_node.at(i, 0) > time_pred){
     if(i == 0)
      temp1 = 1;
     else
      temp1 = leaf_node.at(i-1, 1);
     break;
    } else if (leaf_node.at(i, 0) == time_pred){
     temp1 = leaf_node.at(i, 1);
     break;
    }
   }

  } else {

   // go here if prediction horizon > max time in current leaf.
   temp1 = leaf_node.at(leaf_node.n_rows - 1, 1);

  }

  // running mean: mean_k = mean_{k-1} + (new val - old val) / k
  // compute new val - old val
  // be careful, every oob row has a different denom!
  temp2 = temp1 - surv_pvec[rows_oobag[*iit]];
  surv_pvec[rows_oobag[*iit]] += temp2 / denom_pred[rows_oobag[*iit]];
  ++iit;

  if(iit < iit_vals.end()){

   while(person_leaf == leaf_pred(*iit)){

    temp2 = temp1 - surv_pvec[rows_oobag[*iit]];
    surv_pvec[rows_oobag[*iit]] += temp2 / denom_pred[rows_oobag[*iit]];

    ++iit;

    if (iit == iit_vals.end()) break;

   }

  }

 } while (iit < iit_vals.end());

 // if(verbose > 0){
 //  Rcout << "surv_pvec:" << std::endl << surv_pvec.t() << std::endl;
 // }

}

double oobag_c_harrell(){

 vec time = y_input.unsafe_col(0);
 vec status = y_input.unsafe_col(1);
 iit_vals = find(status == 1);

 k = y_input.n_rows;

 double total=0, concordant=0;

 for (iit = iit_vals.begin(); iit < iit_vals.end(); ++iit) {

  for(j = *iit + 1; j < k; ++j){

   if (time[j] > time[*iit]) { // ties not counted

    total++;

    // since surv_pvec is survival prob (not risk),
    // flip the less than to a greater than.
    if (surv_pvec[j] > surv_pvec[*iit]){

     concordant++;

    } else if (surv_pvec[j] == surv_pvec[*iit]){

     concordant+= 0.5;

    }

   }

  }

 }

 return(concordant / total);

}

// [[Rcpp::export]]
double oobag_c_harrell_testthat(NumericMatrix y_mat,
                                NumericVector s_vec) {

 y_input = mat(y_mat.begin(), y_mat.nrow(), y_mat.ncol(), false);
 surv_pvec = vec(s_vec.begin(), s_vec.length(), false);
 return(oobag_c_harrell());

}

void new_pred_surv_multi_mean(){

 // allocate memory for output
 // surv_pvec.zeros(x_pred.n_rows);

 surv_pvec.set_size(times_pred.size());
 iit_vals = sort_index(leaf_pred, "ascend");
 iit = iit_vals.begin();

 do {

  person_leaf = leaf_pred(*iit);

  for(i = 0; i < leaf_indices.n_rows; i++){
   if(leaf_indices.at(i, 0) == person_leaf){
    break;
   }
  }

  leaf_node = leaf_nodes.rows(leaf_indices(i, 1),
                              leaf_indices(i, 2));

  // if(verbose > 1){
  //  Rcout << "leaf_node:" << std::endl << leaf_node << std::endl;
  // }

  i = 0;

  for(j = 0; j < times_pred.size(); j++){

   time_pred = times_pred(j);

   if(time_pred < leaf_node(leaf_node.n_rows - 1, 0)){

    for(; i < leaf_node.n_rows; i++){

     if (leaf_node(i, 0) > time_pred){

      if(i == 0)
       temp1 = 1;
      else
       temp1 = leaf_node(i-1, 1);

      break;

     } else if (leaf_node(i, 0) == time_pred){

      temp1 = leaf_node(i, 1);
      break;

     }

    }

   } else {

    // go here if prediction horizon > max time in current leaf.
    temp1 = leaf_node(leaf_node.n_rows - 1, 1);

   }

   surv_pvec(j) = temp1;

  }



  surv_pmat.row(*iit) += surv_pvec.t();
  ++iit;

  if(iit < iit_vals.end()){

   while(person_leaf == leaf_pred(*iit)){

    surv_pmat.row(*iit) += surv_pvec.t();
    ++iit;

    if (iit == iit_vals.end()) break;

   }

  }

 } while (iit < iit_vals.end());

}

// void new_pred_surv_multi_median(){
//
//  // allocate memory for output
//  // surv_pvec.zeros(x_pred.n_rows);
//
//  surv_pvec.set_size(times_pred.size());
//  iit_vals = sort_index(leaf_pred, "ascend");
//  iit = iit_vals.begin();
//
//  do {
//
//   person_leaf = leaf_pred(*iit);
//
//   for(i = 0; i < leaf_indices.n_rows; i++){
//    if(leaf_indices.at(i, 0) == person_leaf){
//     break;
//    }
//   }
//
//   leaf_node = leaf_nodes.rows(leaf_indices(i, 1),
//                               leaf_indices(i, 2));
//
//   if(verbose > 1){
//    Rcout << "leaf_node:" << std::endl << leaf_node << std::endl;
//   }
//
//   i = 0;
//
//   for(j = 0; j < times_pred.size(); j++){
//
//    time_pred = times_pred(j);
//
//    if(time_pred < leaf_node(leaf_node.n_rows - 1, 0)){
//
//     for(; i < leaf_node.n_rows; i++){
//
//      if (leaf_node(i, 0) > time_pred){
//
//       if(i == 0)
//        temp1 = 1;
//       else
//        temp1 = leaf_node(i-1, 1);
//
//       break;
//
//      } else if (leaf_node(i, 0) == time_pred){
//
//       temp1 = leaf_node(i, 1);
//       break;
//
//      }
//
//     }
//
//    } else {
//
//     // go here if prediction horizon > max time in current leaf.
//     temp1 = leaf_node(leaf_node.n_rows - 1, 1);
//
//    }
//
//    surv_pvec(j) = temp1;
//
//   }
//
//
//   for(i = 0; i < times_pred.size(); i++){
//    surv_pcube(tree, *iit, i) = surv_pvec[i];
//   }
//
//   ++iit;
//
//   if(iit < iit_vals.end()){
//
//    while(person_leaf == leaf_pred(*iit)){
//
//     for(i = 0; i < times_pred.size(); i++){
//      surv_pcube(tree, *iit, i) = surv_pvec[i];
//     }
//
//     ++iit;
//
//     if (iit == iit_vals.end()) break;
//
//    }
//
//   }
//
//  } while (iit < iit_vals.end());
//
// }


// naming pattern
// new for new data, oob for training
// pred for predict
// surv for survival
// uni for one time, multi for many times
// mean for mean prediction, median for median prediction
void new_pred_surv_uni_mean(){

 iit_vals = sort_index(leaf_pred, "ascend");
 iit = iit_vals.begin();

 do {

  person_leaf = leaf_pred(*iit);

  for(i = 0; i < leaf_indices.n_rows; i++){
   if(leaf_indices.at(i, 0) == person_leaf){
    break;
   }
  }

  leaf_node = leaf_nodes.rows(leaf_indices(i, 1),
                              leaf_indices(i, 2));

  // if(verbose > 1){
  //  Rcout << "leaf_node:" << std::endl << leaf_node << std::endl;
  // }

  i = 0;

  if(time_pred < leaf_node(leaf_node.n_rows - 1, 0)){

   for(; i < leaf_node.n_rows; i++){
    if (leaf_node(i, 0) > time_pred){

     if(i == 0){

      temp1 = 1;

     } else {

      temp1 = leaf_node(i - 1, 1);

      // experimental - does not seem to help!
      // weighted average of surv est from before and after time of pred
      // temp2 = leaf_node(i, 0) - leaf_node(i-1, 0);
      //
      // temp1 = leaf_node(i, 1) * (time_pred - leaf_node(i-1,0)) / temp2 +
      //  leaf_node(i-1, 1) * (leaf_node(i,0) - time_pred) / temp2;

     }

     break;

    } else if (leaf_node(i, 0) == time_pred){
     temp1 = leaf_node(i, 1);
     break;
    }
   }

  } else if (time_pred == leaf_node(leaf_node.n_rows - 1, 0)){

   temp1 = leaf_node(leaf_node.n_rows - 1, 1);

  } else {

   // go here if prediction horizon > max time in current leaf.
   temp1 = leaf_node(leaf_node.n_rows - 1, 1);

   // --- EXPERIMENTAL ADD-ON --- //
   // if you are predicting beyond the max time in a node,
   // then determine how much further out you are and assume
   // the survival probability decays at the same rate.

   // temp2 = (1.0 - temp1) *
   //  (time_pred - leaf_node(leaf_node.n_rows - 1, 0)) / time_pred;
   //
   // temp1 = temp1 * (1.0-temp2);

  }

  // cumulative hazard prediction

  surv_pvec(*iit) += temp1;
  ++iit;

  if(iit < iit_vals.end()){

   while(person_leaf == leaf_pred(*iit)){

    surv_pvec(*iit) += temp1;
    ++iit;

    if (iit == iit_vals.end()) break;

   }

  }

 } while (iit < iit_vals.end());

 // if(verbose > 1){
 //  Rcout << "pred_surv:" << std::endl << surv_pvec.t() << std::endl;
 // }

}

// void new_pred_surv_uni_median(){
//
//  iit_vals = sort_index(leaf_pred, "ascend");
//  iit = iit_vals.begin();
//
//  do {
//
//   person_leaf = leaf_pred(*iit);
//
//   for(i = 0; i < leaf_indices.n_rows; i++){
//    if(leaf_indices.at(i, 0) == person_leaf){
//     break;
//    }
//   }
//
//   leaf_node = leaf_nodes.rows(leaf_indices(i, 1),
//                               leaf_indices(i, 2));
//
//   if(verbose > 1){
//    Rcout << "leaf_node:" << std::endl << leaf_node << std::endl;
//   }
//
//   i = 0;
//
//   if(time_pred < leaf_node(leaf_node.n_rows - 1, 0)){
//
//    for(; i < leaf_node.n_rows; i++){
//     if (leaf_node(i, 0) > time_pred){
//      if(i == 0){
//       temp1 = 1;
//      } else {
//       temp2 = leaf_node(i, 0) - leaf_node(i-1, 0);
//
//       temp1 = leaf_node(i, 1) * (time_pred - leaf_node(i-1,0)) / temp2 +
//        leaf_node(i-1, 1) * (leaf_node(i,0) - time_pred) / temp2;
//      }
//      break;
//     } else if (leaf_node(i, 0) == time_pred){
//      temp1 = leaf_node(i, 1);
//      break;
//     }
//    }
//
//   } else if (time_pred == leaf_node(leaf_node.n_rows - 1, 0)){
//
//    temp1 = leaf_node(leaf_node.n_rows - 1, 1);
//
//   } else {
//
//    // go here if prediction horizon > max time in current leaf.
//    temp1 = leaf_node(leaf_node.n_rows - 1, 1);
//
//    // --- EXPERIMENTAL ADD-ON --- //
//    // if you are predicting beyond the max time in a node,
//    // then determine how much further out you are and assume
//    // the survival probability decays at the same rate.
//
//    temp2 = (1.0 - temp1) *
//     (time_pred - leaf_node(leaf_node.n_rows - 1, 0)) / time_pred;
//
//    temp1 = temp1 * (1.0-temp2);
//
//   }
//
//   surv_pmat(tree, *iit) = temp1;
//   ++iit;
//
//   if(iit < iit_vals.end()){
//
//    while(person_leaf == leaf_pred(*iit)){
//
//     surv_pmat(tree, *iit) = temp1;
//     ++iit;
//
//     if (iit == iit_vals.end()) break;
//
//    }
//
//   }
//
//  } while (iit < iit_vals.end());
//
//
// }


// ----------------------------------------------------------------------------
// --------------------------- ostree functions -------------------------------
// ----------------------------------------------------------------------------

void ostree_size_buffer(){

 // if(verbose > 1){
 //  Rcout << "---------- buffering outputs ----------" << std::endl;
 //  Rcout << "betas before:  " << std::endl << betas.t() << std::endl;
 // }

 betas.insert_cols(betas.n_cols, 10);
 // x_mean.insert_cols(x_mean.n_cols, 10);
 col_indices.insert_cols(col_indices.n_cols, 10);
 children_left.insert_rows(children_left.size(), 10);
 cutpoints.insert_rows(cutpoints.size(), 10);

 // if(verbose > 1){
 //  Rcout << "betas after:  " << std::endl << betas.t() << std::endl;
 //  Rcout << "---------------------------------------";
 //  Rcout << std::endl << std::endl;
 // }


}

void ostree_mem_xfer(){

 // no data copied according to tracemem.
 // not including boot rows or x_mean (don't always need them)

 NumericMatrix leaf_nodes_      = ostree["leaf_nodes"];
 NumericMatrix betas_           = ostree["betas"];
 NumericVector cutpoints_       = ostree["cut_points"];
 IntegerMatrix col_indices_     = ostree["col_indices"];
 IntegerMatrix leaf_indices_    = ostree["leaf_node_index"];
 IntegerVector children_left_   = ostree["children_left"];

 leaf_nodes = mat(leaf_nodes_.begin(),
                  leaf_nodes_.nrow(),
                  leaf_nodes_.ncol(),
                  false);

 betas = mat(betas_.begin(),
             betas_.nrow(),
             betas_.ncol(),
             false);

 cutpoints = vec(cutpoints_.begin(), cutpoints_.length(), false);

 col_indices = conv_to<umat>::from(
  imat(col_indices_.begin(),
       col_indices_.nrow(),
       col_indices_.ncol(),
       false)
 );

 leaf_indices = conv_to<umat>::from(
  imat(leaf_indices_.begin(),
       leaf_indices_.nrow(),
       leaf_indices_.ncol(),
       false)
 );

 children_left = conv_to<uvec>::from(
  ivec(children_left_.begin(),
       children_left_.length(),
       false)
 );

}

void ostree_pred_leaf(){

 // reset values
 // this is needed for pred_leaf since every obs gets a new leaf in
 // the next tree, but it isn't needed for pred_surv because survival
 // probs get aggregated over all the trees.
 leaf_pred.fill(0);

 for(i = 0; i < betas.n_cols; i++){

  if(children_left[i] != 0){

   if(i == 0){
    obs_in_node = regspace<uvec>(0, 1, leaf_pred.size()-1);
   } else {
    obs_in_node = find(leaf_pred == i);
   }


   if(obs_in_node.size() > 0){

    // Fastest sub-matrix multiplication i can think of.
    // Matrix product = linear combination of columns
    // (this is faster b/c armadillo is great at making
    //  pointers to the columns of an arma mat)
    // XB.zeros(obs_in_node.size());
    //
    // uvec col_indices_i = col_indices.unsafe_col(i);
    //
    // j = 0;
    //
    // jit = col_indices_i.begin();
    //
    // for(; jit < col_indices_i.end(); ++jit, ++j){
    //
    //  vec x_j = x_pred.unsafe_col(*jit);
    //
    //  XB += x_j(obs_in_node) * betas.at(j, i);
    //
    // }

    // this is slower but more clear matrix multiplication
    XB = x_pred(obs_in_node, col_indices.col(i)) * betas.col(i);

    jit = obs_in_node.begin();

    for(j = 0; j < XB.size(); ++j, ++jit){

     if(XB[j] <= cutpoints[i]) {

      leaf_pred[*jit] = children_left[i];

     } else {

      leaf_pred[*jit] = children_left[i]+1;

     }

    }

    // if(verbose > 0){
    //
    //  uvec in_left = find(leaf_pred == children_left(i));
    //  uvec in_right = find(leaf_pred == children_left(i)+1);
    //
    //  Rcout << "N to node_" << children_left(i) << ": ";
    //  Rcout << in_left.size() << "; ";
    //  Rcout << "N to node_" << children_left(i)+1 << ": ";
    //  Rcout << in_right.size() << std::endl;
    //
    // }

   }

  }

 }



}

// [[Rcpp::export]]
arma::uvec ostree_pred_leaf_testthat(List& tree,
                                     NumericMatrix& x_pred_){


 x_pred = mat(x_pred_.begin(),
               x_pred_.nrow(),
               x_pred_.ncol(),
               false);

 leaf_pred.set_size(x_pred.n_rows);

 ostree = tree;
 ostree_mem_xfer();
 ostree_pred_leaf();

 return(leaf_pred);

}

List ostree_fit(Function f_beta){

 betas.fill(0);
 // x_mean.fill(0);
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

 // ----------------------
 // ---- main do loop ----
 // ----------------------



 do {

  nodes_to_grow_next.set_size(0);

  // if(verbose > 0){
  //
  //  Rcout << "----------- nodes to grow -----------" << std::endl;
  //  Rcout << "nodes: "<< nodes_to_grow.t()           << std::endl;
  //  Rcout << "-------------------------------------" << std::endl <<
  //   std::endl << std::endl;
  //
  //
  // }

  for(node = nodes_to_grow.begin(); node != nodes_to_grow.end(); ++node){

   if(nodes_to_grow[0] == 0){

    // when growing the first node, there is no need to find
    // which rows are in the node.
    rows_node = linspace<uvec>(0,
                               x_inbag.n_rows-1,
                               x_inbag.n_rows);

   } else {

    // identify which rows are in the current node.
    rows_node = find(node_assignments == *node);

   }

   y_node = y_inbag.rows(rows_node);
   w_node = w_inbag(rows_node);

   // if(verbose > 0){
   //
   //  n_risk = sum(w_node);
   //  n_events = sum(y_node.col(1) % w_node);
   //  Rcout << "-------- Growing node " << *node << " --------" << std::endl;
   //  Rcout << "No. of observations in node: " << n_risk        << std::endl;
   //  Rcout << "No. of events in node:       " << n_events      << std::endl;
   //  Rcout << "No. of rows in node:         " << w_node.size() << std::endl;
   //  Rcout << "--------------------------------"               << std::endl;
   //  Rcout << std::endl << std::endl;
   //
   // }

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

   n_cols_to_sample = sum(cols_to_sample_01);

   if(n_cols_to_sample > 1){

    n_events_total = sum(y_node.col(1) % w_node);

    if(n_cols_to_sample < mtry){

     mtry_int = n_cols_to_sample;

     // if(verbose > 0){
     //  Rcout << " ---- >=1 constant column in node rows ----" << std::endl;
     //  Rcout << "mtry reduced to " << mtry_temp << " from " << mtry;
     //  Rcout << std::endl;
     //  Rcout << "-------------------------------------------" << std::endl;
     //  Rcout << std::endl << std::endl;
     // }

    }

    if (type_beta == 'C'){

     // make sure there are at least 3 event per predictor variable.
     // (if using CPH)
     while(n_events_total / mtry_int < 3 && mtry_int > 1){
      --mtry_int;
     }

    }


    n_cols_to_sample = mtry_int;

    // if(verbose > 0){
    //  Rcout << "n_events: " << n_events_total << std::endl;
    //  Rcout << "mtry: " << mtry_int << std::endl;
    //  Rcout << "n_events per column: " << n_events_total/mtry_int << std::endl;
    // }

    if(mtry_int > 1){

     cols_to_sample = find(cols_to_sample_01);

     // re-try hinge point
     n_retry = 0;
     cutpoint = R_PosInf;

     while(n_retry <= max_retry){

      // if(n_retry > 0) Rcout << "trying again!" << std::endl;

      cols_node = Rcpp::RcppArmadillo::sample(cols_to_sample,
                                              mtry_int,
                                              false);

      x_node = x_inbag(rows_node, cols_node);

      // here is where n_vars gets updated to match the current node
      // originally it matched the number of variables in the input x.

      n_vars = x_node.n_cols;

      if(cph_do_scale){
       x_node_scale();
      }

      // if(verbose > 0){
      //
      //  uword temp_uword_1 = min(uvec {x_node.n_rows, 5});
      //  Rcout << "x node scaled: " << std::endl;
      //  Rcout << x_node.submat(0, 0, temp_uword_1-1, x_node.n_cols-1);
      //  Rcout << std::endl;
      //
      // }

      switch(type_beta) {

      case 'C' :

       beta_fit = newtraph_cph();

       if(cph_do_scale){
        for(i = 0; i < x_transforms.n_rows; i++){
         x_node.col(i) /= x_transforms(i,1);
         x_node.col(i) += x_transforms(i,0);
        }

       }

       break;

      case 'N' :

       xx = wrap(x_node);
       yy = wrap(y_node);
       ww = wrap(w_node);
       colnames(yy) = yy_names;

       beta_placeholder = f_beta(xx, yy, ww,
                                 net_alpha,
                                 net_df_target);

       beta_fit = mat(beta_placeholder.begin(),
                      beta_placeholder.nrow(),
                      beta_placeholder.ncol(),
                      false);

       break;

      case 'U' :

       xx = wrap(x_node);
       yy = wrap(y_node);
       ww = wrap(w_node);
       colnames(yy) = yy_names;

       beta_placeholder = f_beta(xx, yy, ww);

       beta_fit = mat(beta_placeholder.begin(),
                      beta_placeholder.nrow(),
                      beta_placeholder.ncol(),
                      false);

       break;

      }


      if(any(beta_fit)){

       // if(verbose > 0){
       //
       //  uword temp_uword_1 = min(uvec {x_node.n_rows, 5});
       //  Rcout << "x node unscaled: " << std::endl;
       //  Rcout << x_node.submat(0, 0, temp_uword_1-1, x_node.n_cols-1);
       //  Rcout << std::endl;
       //
       // }

       XB = x_node * beta_fit;
       cutpoint = lrt_multi();

      }

      if(!std::isinf(cutpoint)) break;
      n_retry++;

     }

    }

   }

   if(!std::isinf(cutpoint)){

    // make new nodes if a valid cutpoint was found
    nn_left   = nodes_max_true + 1;
    nodes_max_true = nodes_max_true + 2;


    // if(verbose > 0){
    //
    //  Rcout << "-------- New nodes created --------" << std::endl;
    //  Rcout << "Left node: "  << nn_left             << std::endl;
    //  Rcout << "Right node: " << nodes_max_true      << std::endl;
    //  Rcout << "-----------------------------------" << std::endl <<
    //   std::endl << std::endl;
    //
    // }

    n_events_left = n_events_total - n_events_right;

    // if(verbose > 0){
    //  Rcout << "n_events_left: " << n_events_left << std::endl;
    //  Rcout << "n_risk_left: " << n_risk_left << std::endl;
    //  Rcout << "n_events_right: " << n_events_right << std::endl;
    //  Rcout << "n_risk_right: " << n_risk_right << std::endl;
    // }

    i=0;

    for(iit = rows_node.begin(); iit != rows_node.end(); ++iit, ++i){

     node_assignments[*iit] = nn_left + group[i];

    }

    if(n_events_left >= 2*leaf_min_events &&
       n_risk_left   >= 2*leaf_min_obs &&
       n_events_left >=   split_min_events &&
       n_risk_left   >=   split_min_obs){

     nodes_to_grow_next = join_cols(nodes_to_grow_next,
                                    uvec{nn_left});

    } else {

     rows_leaf = find(group==0);
     leaf_indices(leaf_node_index_counter, 0) = nn_left;
     leaf_kaplan(y_node.rows(rows_leaf), w_node(rows_leaf));

     // if(verbose > 0){
     //  Rcout << "-------- creating a new leaf --------" << std::endl;
     //  Rcout << "name: node_" << nn_left                << std::endl;
     //  Rcout << "n_obs:    "  << sum(w_node(rows_leaf));
     //  Rcout << std::endl;
     //  Rcout << "n_events: ";
     //  vec_temp = y_node.col(1);
     //  Rcout << sum(w_node(rows_leaf) % vec_temp(rows_leaf));
     //  Rcout << std::endl;
     //  Rcout << "------------------------------------";
     //  Rcout << std::endl << std::endl << std::endl;
     // }

    }

    if(n_events_right >= 2*leaf_min_events &&
       n_risk_right   >= 2*leaf_min_obs &&
       n_events_right >=   split_min_events &&
       n_risk_right   >=   split_min_obs){

     nodes_to_grow_next = join_cols(nodes_to_grow_next,
                                    uvec{nodes_max_true});

    } else {

     rows_leaf = find(group==1);
     leaf_indices(leaf_node_index_counter, 0) = nodes_max_true;
     leaf_kaplan(y_node.rows(rows_leaf), w_node(rows_leaf));

     // if(verbose > 0){
     //  Rcout << "-------- creating a new leaf --------" << std::endl;
     //  Rcout << "name: node_" << nodes_max_true               << std::endl;
     //  Rcout << "n_obs:    "  << sum(w_node(rows_leaf));
     //  Rcout << std::endl;
     //  Rcout << "n_events: ";
     //  vec_temp = y_node.col(1);
     //  Rcout << sum(w_node(rows_leaf) % vec_temp(rows_leaf));
     //  Rcout << std::endl;
     //  Rcout << "------------------------------------";
     //  Rcout << std::endl << std::endl << std::endl;
     // }

    }

    if(nodes_max_true >= betas.n_cols) ostree_size_buffer();

    for(i = 0; i < n_cols_to_sample; i++){
     betas.at(i, *node) = beta_fit[i];
     // x_mean.at(i, *node) = x_transforms(i, 0);
     col_indices.at(i, *node) = cols_node[i];
    }

    children_left[*node] = nn_left;
    cutpoints[*node] = cutpoint;

   } else {

    // make a leaf node if a valid cutpoint could not be found
    leaf_indices(leaf_node_index_counter, 0) = *node;
    leaf_kaplan(y_node, w_node);

    // if(verbose > 0){
    //  Rcout << "-------- creating a new leaf --------" << std::endl;
    //  Rcout << "name: node_" << *node                  << std::endl;
    //  Rcout << "n_obs:    "  << sum(w_node)      << std::endl;
    //  Rcout << "n_events: "  << sum(w_node % y_node.col(1));
    //  Rcout                                            << std::endl;
    //  Rcout << "Couldn't find a cutpoint??"            << std::endl;
    //  Rcout << "------------------------------------"  << std::endl;
    //  Rcout << std::endl << std::endl;
    // }

   }

  }

  nodes_to_grow = nodes_to_grow_next;

 } while (nodes_to_grow.size() > 0);

 return(
  List::create(

   _["leaf_nodes"] = leaf_nodes.rows(span(0, leaf_node_counter-1)),

   _["leaf_node_index"] = conv_to<imat>::from(
    leaf_indices.rows(span(0, leaf_node_index_counter-1))
   ),

   _["betas"] = betas.cols(span(0, nodes_max_true)),

   // _["x_mean"] = x_mean.cols(span(0, nodes_max_true)),

   _["col_indices"] = conv_to<imat>::from(
    col_indices.cols(span(0, nodes_max_true))
   ),

   _["cut_points"] = cutpoints(span(0, nodes_max_true)),

   _["children_left"] = conv_to<ivec>::from(
    children_left(span(0, nodes_max_true))
   ),

   _["rows_oobag"] = conv_to<ivec>::from(rows_oobag)

  )
 );


}

// ----------------------------------------------------------------------------
// ---------------------------- orsf functions --------------------------------
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
List orsf_fit(NumericMatrix& x,
              NumericMatrix& y,
              const int&     n_tree,
              const int&     n_split_,
              const int&     mtry_,
              const double&  leaf_min_events_,
              const double&  leaf_min_obs_,
              const double&  split_min_events_,
              const double&  split_min_obs_,
              const double&  split_min_stat_,
              const int&     cph_method_,
              const double&  cph_eps_,
              const int&     cph_iter_max_,
              const bool&    cph_do_scale_,
              const double&  net_alpha_,
              const int&     net_df_target_,
              const bool&    oobag_pred_,
              const double&  oobag_pred_horizon_,
              const int&     oobag_eval_every_,
              const bool&    oobag_importance_,
              const char&    oobag_importance_type_,
              IntegerVector& tree_seeds,
              const int&     max_retry_,
              Function       f_beta,
              const char&    type_beta_,
              Function       f_oobag_eval,
              const char&    type_oobag_eval_){


 // convert inputs into arma objects
 x_input = mat(x.begin(), x.nrow(), x.ncol(), false);
 y_input = mat(y.begin(), y.nrow(), y.ncol(), false);

 // these change later in ostree_fit()
 n_rows = x_input.n_rows;
 n_vars = x_input.n_cols;

 // initialize the variable importance (vi) vectors
 vi_pval_numer.zeros(n_vars);
 vi_pval_denom.zeros(n_vars);

 // if(verbose > 0){
 //  Rcout << "------------ dimensions ------------"  << std::endl;
 //  Rcout << "N obs total: "     << n_rows           << std::endl;
 //  Rcout << "N columns total: " << n_vars           << std::endl;
 //  Rcout << "------------------------------------";
 //  Rcout << std::endl << std::endl << std::endl;
 // }

 n_split               = n_split_;
 mtry                  = mtry_;
 leaf_min_events       = leaf_min_events_;
 leaf_min_obs          = leaf_min_obs_;
 split_min_events      = split_min_events_;
 split_min_obs         = split_min_obs_;
 split_min_stat        = split_min_stat_;
 cph_method            = cph_method_;
 cph_eps               = cph_eps_;
 cph_iter_max          = cph_iter_max_;
 cph_do_scale          = cph_do_scale_;
 net_alpha             = net_alpha_;
 net_df_target         = net_df_target_;
 oobag_pred            = oobag_pred_;
 oobag_eval_every      = oobag_eval_every_;
 oobag_eval_counter    = 0;
 oobag_importance      = oobag_importance_;
 oobag_importance_type = oobag_importance_type_;
 use_tree_seed         = tree_seeds.length() > 0;
 max_retry             = max_retry_;
 type_beta             = type_beta_;
 type_oobag_eval       = type_oobag_eval_;
 temp1                 = 1.0 / n_rows;

 if(cph_iter_max > 1) cph_do_scale = true;

 if((type_beta == 'N') || (type_beta == 'U')) cph_do_scale = false;

 if(oobag_pred){

  time_pred = oobag_pred_horizon_;

  if(time_pred == 0) time_pred = median(y_input.col(0));

  eval_oobag.set_size(std::floor(n_tree / oobag_eval_every));

 } else {

  eval_oobag.set_size(0);

 }

 // if(verbose > 0){
 //  Rcout << "------------ input variables ------------" << std::endl;
 //  Rcout << "n_split: "         << n_split              << std::endl;
 //  Rcout << "mtry: "            << mtry                 << std::endl;
 //  Rcout << "leaf_min_events: " << leaf_min_events      << std::endl;
 //  Rcout << "leaf_min_obs: "    << leaf_min_obs         << std::endl;
 //  Rcout << "cph_method: "      << cph_method           << std::endl;
 //  Rcout << "cph_eps: "         << cph_eps              << std::endl;
 //  Rcout << "cph_iter_max: "    << cph_iter_max         << std::endl;
 //  Rcout << "-----------------------------------------" << std::endl;
 //  Rcout << std::endl << std::endl;
 // }

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

  surv_pvec.zeros(n_rows);
  denom_pred.zeros(n_rows);

 } else {

  surv_pvec.set_size(0);
  denom_pred.set_size(0);

 }

 // guessing the number of nodes needed to grow a tree
 nodes_max_guess = std::ceil(0.5 * n_rows / leaf_min_events);

 betas.zeros(mtry, nodes_max_guess);
 // x_mean.zeros(mtry, nodes_max_guess);
 col_indices.zeros(mtry, nodes_max_guess);
 cutpoints.zeros(nodes_max_guess);
 children_left.zeros(nodes_max_guess);
 leaf_indices.zeros(nodes_max_guess, 3);

 // some great variable names here
 List forest(n_tree);

 for(tree = 0; tree < n_tree; ){

  // Abort the routine if user has pressed Ctrl + C or Escape in R.
  Rcpp::checkUserInterrupt();

  // --------------------------------------------
  // ---- initialize parameters to grow tree ----
  // --------------------------------------------

  // rows_inbag = find(w_inbag != 0);

  if(use_tree_seed) set_seed_r(tree_seeds[tree]);

  w_input    = as<vec>(sample(s, n_rows, true, probs));
  rows_oobag = find(w_input == 0);
  rows_inbag = regspace<uvec>(0, n_rows-1);
  rows_inbag = std_setdiff(rows_inbag, rows_oobag);
  w_inbag    = w_input(rows_inbag);

  // if(verbose > 0){
  //
  //  Rcout << "------------ boot weights ------------" << std::endl;
  //  Rcout << "pr(inbag): " << 1-pow(1-temp1,n_rows)   << std::endl;
  //  Rcout << "total: "     << sum(w_inbag)      << std::endl;
  //  Rcout << "N > 0: "     << rows_inbag.size()       << std::endl;
  //  Rcout << "--------------------------------------" <<
  //   std::endl << std::endl << std::endl;
  //
  // }

  x_inbag = x_input.rows(rows_inbag);
  y_inbag = y_input.rows(rows_inbag);

  if(oobag_pred){
   x_pred = x_input.rows(rows_oobag);
   leaf_pred.set_size(rows_oobag.size());
  }

  // if(verbose > 0){
  //
  //  uword temp_uword_1, temp_uword_2;
  //
  //  if(x_inbag.n_rows < 5)
  //   temp_uword_1 = x_inbag.n_rows-1;
  //  else
  //   temp_uword_1 = 5;
  //
  //  if(x_inbag.n_cols < 5)
  //   temp_uword_2 = x_inbag.n_cols-1;
  //  else
  //   temp_uword_2 = 4;
  //
  //  Rcout << "x inbag: " << std::endl <<
  //   x_inbag.submat(0, 0,
  //                  temp_uword_1,
  //                  temp_uword_2) << std::endl;
  //
  // }

  forest[tree] = ostree_fit(f_beta);

  // add 1 to tree here instead of end of loop
  // (more convenient to compute tree % oobag_eval_every)
  tree++;

  if(oobag_pred){

   denom_pred(rows_oobag) += 1;
   ostree_pred_leaf();
   oobag_pred_surv_uni();

   if(tree % oobag_eval_every == 0){

    switch(type_oobag_eval) {

    // H stands for Harrell's C-statistic
    case 'H' :

     eval_oobag[oobag_eval_counter] = oobag_c_harrell();
     oobag_eval_counter++;

     break;

    // U stands for a user-supplied function
    case 'U' :

     ww = wrap(surv_pvec);

     eval_oobag[oobag_eval_counter] = as<double>(
      f_oobag_eval(y, ww)
     );

     oobag_eval_counter++;

     break;

    }


   }

  }

 }



 vec vimp(x_input.n_cols);

 // ANOVA importance
 if(oobag_importance_type == 'A') vimp = vi_pval_numer / vi_pval_denom;

 // if we are computing variable importance, surv_pvec is about
 // to get modified, and we don't want to return the modified
 // version of surv_pvec.
 // So make a deep copy if oobag_importance is true.
 // Make a shallow copy if oobag_importance is false
 surv_pvec_output = vec(surv_pvec.begin(),
                        surv_pvec.size(),
                        oobag_importance);

 if(oobag_importance && n_tree > 0){

  uvec betas_to_flip;
  oobag_eval_counter--;

  for(uword variable = 0; variable < x_input.n_cols; ++variable){

   surv_pvec.fill(0);
   denom_pred.fill(0);

   for(tree = 0; tree < n_tree; ++tree){

    ostree = forest[tree];

    IntegerMatrix rows_oobag_ = ostree["rows_oobag"];

    rows_oobag = conv_to<uvec>::from(
     ivec(rows_oobag_.begin(),
          rows_oobag_.length(),
          false)
    );

    x_pred = x_input.rows(rows_oobag);

    if(oobag_importance_type == 'P'){
     x_pred.col(variable) = shuffle(x_pred.col(variable));
    }

    ostree_mem_xfer();


    if(oobag_importance_type == 'N'){
     betas_to_flip = find(col_indices == variable);
     betas.elem( betas_to_flip ) *= (-1);
    }

    denom_pred(rows_oobag) += 1;

    leaf_pred.set_size(rows_oobag.size());

    ostree_pred_leaf();

    oobag_pred_surv_uni();

    if(oobag_importance_type == 'N'){
     betas.elem( betas_to_flip ) *= (-1);
    }

   }

   switch(type_oobag_eval) {

   // H stands for Harrell's C-statistic
   case 'H' :

    vimp(variable) = eval_oobag[oobag_eval_counter] - oobag_c_harrell();

    break;

    // U stands for a user-supplied function
   case 'U' :

    ww = wrap(surv_pvec);

    vimp(variable) =
     eval_oobag[oobag_eval_counter] - as<double>(f_oobag_eval(y, ww));


    break;

   }

  }

 }

 return(
  List::create(
   _["forest"] = forest,
   _["surv_oobag"] = surv_pvec_output,
   _["pred_horizon"] = time_pred,
   _["eval_oobag"] = List::create(_["stat_values"] = eval_oobag,
                                  _["stat_type"]   = type_oobag_eval),
   _["importance"] = vimp
  )
 );


}

// [[Rcpp::export]]
arma::vec orsf_oob_negate_vi(NumericMatrix& x,
                             NumericMatrix& y,
                             List& forest,
                             const double& last_eval_stat,
                             const double& time_pred_,
                             Function      f_oobag_eval,
                             const char&   type_oobag_eval_){

 x_input = mat(x.begin(), x.nrow(), x.ncol(), false);
 y_input = mat(y.begin(), y.nrow(), y.ncol(), false);

 time_pred = time_pred_;
 type_oobag_eval = type_oobag_eval_;

 vec vimp(x_input.n_cols);

 uvec betas_to_flip;
 uword variable;

 for(variable = 0; variable < x_input.n_cols; ++variable){

  surv_pvec.fill(0);
  denom_pred.fill(0);

  for(tree = 0; tree < forest.length(); ++tree){

   ostree = forest[tree];

   IntegerMatrix rows_oobag_ = ostree["rows_oobag"];

   rows_oobag = conv_to<uvec>::from(
    ivec(rows_oobag_.begin(),
         rows_oobag_.length(),
         false)
   );

   x_pred = x_input.rows(rows_oobag);

   ostree_mem_xfer();

   betas_to_flip = find(col_indices == variable);

   betas.elem( betas_to_flip ) *= (-1);

   denom_pred(rows_oobag) += 1;

   leaf_pred.set_size(rows_oobag.size());

   ostree_pred_leaf();

   oobag_pred_surv_uni();

   betas.elem( betas_to_flip ) *= (-1);

  }

  switch(type_oobag_eval) {

  // H stands for Harrell's C-statistic
  case 'H' :

   vimp(variable) = last_eval_stat - oobag_c_harrell();

   break;

   // U stands for a user-supplied function
  case 'U' :

   ww = wrap(surv_pvec);

   vimp(variable) = last_eval_stat - as<double>(f_oobag_eval(y, ww));

   break;

  }

 }

 return(vimp);

}

// [[Rcpp::export]]
arma::vec orsf_oob_permute_vi(NumericMatrix& x,
                              NumericMatrix& y,
                              List& forest,
                              const double& last_eval_stat,
                              const double& time_pred_,
                              Function      f_oobag_eval,
                              const char&   type_oobag_eval_){

 x_input = mat(x.begin(), x.nrow(), x.ncol(), false);
 y_input = mat(y.begin(), y.nrow(), y.ncol(), false);

 time_pred = time_pred_;
 type_oobag_eval = type_oobag_eval_;

 vec vimp(x_input.n_cols);

 uword variable;

 for(variable = 0; variable < x_input.n_cols; ++variable){

  surv_pvec.fill(0);
  denom_pred.fill(0);

  for(tree = 0; tree < forest.length(); ++tree){

   ostree = forest[tree];

   IntegerMatrix rows_oobag_ = ostree["rows_oobag"];

   rows_oobag = conv_to<uvec>::from(
    ivec(rows_oobag_.begin(),
         rows_oobag_.length(),
         false)
   );

   x_pred = x_input.rows(rows_oobag);

   x_pred.col(variable) = shuffle(x_pred.col(variable));

   ostree_mem_xfer();

   denom_pred(rows_oobag) += 1;

   leaf_pred.set_size(rows_oobag.size());

   ostree_pred_leaf();

   oobag_pred_surv_uni();

   // x_variable = x_variable_original;
   // x_input.col(variable) = x_variable;

  }

  switch(type_oobag_eval) {

  // H stands for Harrell's C-statistic
  case 'H' :

   vimp(variable) = last_eval_stat - oobag_c_harrell();

   break;

   // U stands for a user-supplied function
  case 'U' :

   ww = wrap(surv_pvec);

   vimp(variable) = last_eval_stat - as<double>(f_oobag_eval(y, ww));

   break;

  }

 }

 return(vimp);

}

// [[Rcpp::export]]
arma::mat orsf_pred_uni(List& forest,
                        NumericMatrix& x_new,
                        double time_dbl,
                        bool return_risk){

 x_pred = mat(x_new.begin(), x_new.nrow(), x_new.ncol(), false);
 time_pred = time_dbl;

 // memory for outputs
 leaf_pred.set_size(x_pred.n_rows);
 surv_pvec.zeros(x_pred.n_rows);

 // if(compute_median){
 //
 //  surv_pmat.zeros(forest.length(), x_pred.n_rows);
 //
 //  for(tree = 0; tree < forest.length(); ++tree){
 //   ostree = forest[tree];
 //   ostree_mem_xfer();
 //   ostree_pred_leaf();
 //   new_pred_surv_uni_median();
 //  }
 //
 //  for(i = 0; i < x_pred.n_rows; i++){
 //   surv_pvec[i] = median(surv_pmat.col(i));
 //  }
 //
 // }

  for(tree = 0; tree < forest.length(); ++tree){
   ostree = forest[tree];
   ostree_mem_xfer();
   ostree_pred_leaf();
   new_pred_surv_uni_mean();
  }

  surv_pvec /= tree;

 if(return_risk){
  return(1 - surv_pvec);
 } else {
  return(surv_pvec);
 }

}

// [[Rcpp::export]]
arma::mat orsf_pred_multi(List& forest,
                          NumericMatrix& x_new,
                          NumericVector& time_vec,
                          bool return_risk){

 x_pred = mat(x_new.begin(), x_new.nrow(), x_new.ncol(), false);
 times_pred = vec(time_vec.begin(), time_vec.length(), false);

 // memory for outputs
 // initial values don't matter for leaf_pred,
 // but do matter for surv_pmat
 leaf_pred.set_size(x_pred.n_rows);
 surv_pmat.zeros(x_pred.n_rows, times_pred.size());

 // if(compute_median){
 //
 //  surv_pcube.zeros(forest.length(), x_pred.n_rows, times_pred.size());
 //
 //  for(tree = 0; tree < forest.length(); ++tree){
 //   ostree = forest[tree];
 //   ostree_mem_xfer();
 //   ostree_pred_leaf();
 //   new_pred_surv_multi_median();
 //  }
 //
 //  for(i = 0; i < x_pred.n_rows; i++){
 //   for(j = 0; j < times_pred.size(); j++){
 //    surv_pmat(i, j) = median(surv_pcube.slice(j).col(i));
 //   }
 //  }
 //
 // }


  for(tree = 0; tree < forest.length(); ++tree){
   ostree = forest[tree];
   ostree_mem_xfer();
   ostree_pred_leaf();
   new_pred_surv_multi_mean();
  }

  surv_pmat /= tree;

 if(return_risk){
  return(1 - surv_pmat);
 } else {
  return(surv_pmat);
 }

}

// [[Rcpp::export]]
arma::mat pd_new_smry(List&          forest,
                      NumericMatrix& x_new_,
                      IntegerVector& x_cols_,
                      NumericMatrix& x_vals_,
                      NumericVector& probs_,
                      const double   time_dbl,
                      const bool     return_risk){


 uword pd_i;

 time_pred = time_dbl;

 x_pred = mat(x_new_.begin(), x_new_.nrow(), x_new_.ncol(), false);

 mat x_vals = mat(x_vals_.begin(), x_vals_.nrow(), x_vals_.ncol(), false);

 uvec x_cols = conv_to<uvec>::from(
  ivec(x_cols_.begin(), x_cols_.length(), false)
 );

 vec probs = vec(probs_.begin(), probs_.length(), false);

 mat output_quantiles(probs.size(), x_vals.n_rows);
 mat output_means(1, x_vals.n_rows);

 leaf_pred.set_size(x_pred.n_rows);
 surv_pvec.set_size(x_pred.n_rows);

 for(pd_i = 0; pd_i < x_vals.n_rows; pd_i++){

  j = 0;

  surv_pvec.fill(0);

  for(jit = x_cols.begin(); jit < x_cols.end(); ++jit, ++j){

   x_pred.col(*jit).fill(x_vals(pd_i, j));

  }

  for(tree = 0; tree < forest.length(); ++tree){
   ostree = forest[tree];
   ostree_mem_xfer();
   ostree_pred_leaf();
   new_pred_surv_uni_mean();
  }

  surv_pvec /= tree;

  if(return_risk){ surv_pvec = 1 - surv_pvec; }

  output_means.col(pd_i) = mean(surv_pvec);
  output_quantiles.col(pd_i) = quantile(surv_pvec, probs);


 }

 return(join_vert(output_means, output_quantiles));

}



// [[Rcpp::export]]
arma::mat pd_oob_smry(List&          forest,
                      NumericMatrix& x_new_,
                      IntegerVector& x_cols_,
                      NumericMatrix& x_vals_,
                      NumericVector& probs_,
                      const double   time_dbl,
                      const bool     return_risk){


 uword pd_i;

 time_pred = time_dbl;

 mat x_vals = mat(x_vals_.begin(), x_vals_.nrow(), x_vals_.ncol(), false);

 uvec x_cols = conv_to<uvec>::from(
  ivec(x_cols_.begin(), x_cols_.length(), false)
 );

 vec probs = vec(probs_.begin(), probs_.length(), false);

 mat output_quantiles(probs.size(), x_vals.n_rows);
 mat output_means(1, x_vals.n_rows);

 x_input = mat(x_new_.begin(), x_new_.nrow(), x_new_.ncol(), false);
 denom_pred.set_size(x_input.n_rows);
 surv_pvec.set_size(x_input.n_rows);

 for(pd_i = 0; pd_i < x_vals.n_rows; pd_i++){

  j = 0;
  denom_pred.fill(0);
  surv_pvec.fill(0);

  for(jit = x_cols.begin(); jit < x_cols.end(); ++jit, ++j){

   x_input.col(*jit).fill(x_vals(pd_i, j));

  }

  for(tree = 0; tree < forest.length(); ++tree){

   ostree = forest[tree];

   IntegerMatrix rows_oobag_ = ostree["rows_oobag"];

   rows_oobag = conv_to<uvec>::from(
    ivec(rows_oobag_.begin(),
         rows_oobag_.length(),
         false)
   );

   x_pred = x_input.rows(rows_oobag);
   leaf_pred.set_size(x_pred.n_rows);
   denom_pred(rows_oobag) += 1;

   ostree_mem_xfer();
   ostree_pred_leaf();
   oobag_pred_surv_uni();


  }

  if(return_risk){ surv_pvec = 1 - surv_pvec; }

  output_means.col(pd_i) = mean(surv_pvec);
  output_quantiles.col(pd_i) = quantile(surv_pvec, probs);


 }


 return(join_vert(output_means, output_quantiles));

}

// [[Rcpp::export]]
arma::mat pd_new_ice(List&          forest,
                     NumericMatrix& x_new_,
                     IntegerVector& x_cols_,
                     NumericMatrix& x_vals_,
                     NumericVector& probs_,
                     const double   time_dbl,
                     const bool     return_risk){


 uword pd_i;

 time_pred = time_dbl;

 x_pred = mat(x_new_.begin(), x_new_.nrow(), x_new_.ncol(), false);

 mat x_vals = mat(x_vals_.begin(), x_vals_.nrow(), x_vals_.ncol(), false);

 uvec x_cols = conv_to<uvec>::from(
  ivec(x_cols_.begin(), x_cols_.length(), false)
 );

 vec probs = vec(probs_.begin(), probs_.length(), false);

 mat output_ice(x_vals.n_rows * x_pred.n_rows, 2);
 vec output_ids = output_ice.unsafe_col(0);
 vec output_pds = output_ice.unsafe_col(1);

 uvec pd_rows = regspace<uvec>(0, 1, x_pred.n_rows - 1);

 leaf_pred.set_size(x_pred.n_rows);
 surv_pvec.set_size(x_pred.n_rows);

 for(pd_i = 0; pd_i < x_vals.n_rows; pd_i++){

  j = 0;

  surv_pvec.fill(0);

  for(jit = x_cols.begin(); jit < x_cols.end(); ++jit, ++j){

   x_pred.col(*jit).fill(x_vals(pd_i, j));

  }

  for(tree = 0; tree < forest.length(); ++tree){
   ostree = forest[tree];
   ostree_mem_xfer();
   ostree_pred_leaf();
   new_pred_surv_uni_mean();
  }

  surv_pvec /= tree;

  if(return_risk){ surv_pvec = 1 - surv_pvec; }

  output_ids(pd_rows).fill(pd_i+1);
  output_pds(pd_rows) = surv_pvec;
  pd_rows += x_pred.n_rows;


 }

 return(output_ice);

}

// [[Rcpp::export]]
arma::mat pd_oob_ice(List&          forest,
                     NumericMatrix& x_new_,
                     IntegerVector& x_cols_,
                     NumericMatrix& x_vals_,
                     NumericVector& probs_,
                     const double   time_dbl,
                     const bool     return_risk){


 uword pd_i;

 time_pred = time_dbl;

 mat x_vals = mat(x_vals_.begin(), x_vals_.nrow(), x_vals_.ncol(), false);

 uvec x_cols = conv_to<uvec>::from(
  ivec(x_cols_.begin(), x_cols_.length(), false)
 );

 x_input = mat(x_new_.begin(), x_new_.nrow(), x_new_.ncol(), false);

 mat output_ice(x_vals.n_rows * x_input.n_rows, 2);
 vec output_ids = output_ice.unsafe_col(0);
 vec output_pds = output_ice.unsafe_col(1);

 uvec pd_rows = regspace<uvec>(0, 1, x_input.n_rows - 1);

 denom_pred.set_size(x_input.n_rows);
 surv_pvec.set_size(x_input.n_rows);

 for(pd_i = 0; pd_i < x_vals.n_rows; pd_i++){

  j = 0;
  denom_pred.fill(0);
  surv_pvec.fill(0);

  for(jit = x_cols.begin(); jit < x_cols.end(); ++jit, ++j){

   x_input.col(*jit).fill(x_vals(pd_i, j));

  }

  for(tree = 0; tree < forest.length(); ++tree){

   ostree = forest[tree];

   IntegerMatrix rows_oobag_ = ostree["rows_oobag"];

   rows_oobag = conv_to<uvec>::from(
    ivec(rows_oobag_.begin(),
         rows_oobag_.length(),
         false)
   );

   x_pred = x_input.rows(rows_oobag);
   leaf_pred.set_size(x_pred.n_rows);
   denom_pred(rows_oobag) += 1;

   ostree_mem_xfer();
   ostree_pred_leaf();
   oobag_pred_surv_uni();


  }

  if(return_risk){ surv_pvec = 1 - surv_pvec; }

  output_ids(pd_rows).fill(pd_i+1);
  output_pds(pd_rows) = surv_pvec;
  pd_rows += x_input.n_rows;


 }

 return(output_ice);

}



