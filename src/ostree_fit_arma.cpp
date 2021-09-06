#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat leaf_surv_small(const arma::mat& y,
                          const arma::uvec& weights){

  arma::uword n_dead, n_cens, n_risk, n_risk_sub, person, km_counter;
  arma::colvec time_unique(y.n_rows);

  // use sorted y times to count the number of unique.
  // also define number at risk as the sum of the weights
  arma::uword n_unique = 1; // set at 1 to account for the first time

  // find the first unique event time
  person = 0;
  n_risk = 0;

  while(y(person, 1) == 0){
    n_risk += weights[person];
    person++;
  }

  // now person should correspond to the first event time
  time_unique(0) = y(person,0);  // see above
  double anchor = y(person, 0);

  for( ; person < y.n_rows; person++){


    if(anchor != y(person,0) && y(person,1) == 1){
      time_unique(n_unique) = y(person,0);
      anchor = y(person, 0);
      n_unique++;
    }

    n_risk += weights[person];

  }

  // drop the extra zeros from time_unique
  time_unique = time_unique(arma::span(0, n_unique-1));
  //Rcout << n_risk << std::endl;
  //Rcout << n_unique << std::endl;

  // // reset for next loop
  person = 0;
  km_counter = 0;
  double km = 1.0;
  arma::colvec kmvec(n_unique);

  do {

    anchor = y(person, 0);
    n_dead = 0;
    n_cens = 0;
    n_risk_sub = 0;

    while(y(person, 0) == anchor){

      n_risk_sub += weights[person];

      if(y(person, 1) == 1){
        n_dead += weights[person];
      } else {
        n_cens += weights[person];
      }

      if(person == y.n_rows-1) break;
      person++;

    }

    //Rcout << "n_risk: " << n_risk << std::endl;
    //Rcout << "n_dead: " << n_dead << std::endl;
    //Rcout << "n_risk: " << n_risk << std::endl;

    // only do km if a death was observed

    if(n_dead > 0){

      km = km * (n_risk - n_dead) / n_risk;

      //Rcout << "km: " << km << std::endl;

      kmvec[km_counter] = km;


      km_counter++;

    }

    n_risk -= n_risk_sub;

  } while (km_counter < n_unique);

  //Rcout << time_unique.n_rows << std::endl;
  //Rcout << kmvec.n_rows << std::endl;

  return(arma::join_horiz(time_unique, kmvec));

}

// [[Rcpp::export]]
arma::mat leaf_surv(arma::mat& y,
                    arma::uvec& weights){

  arma::uword n_dead, n_cens, n_risk, n_risk_sub, person, km_counter;
  arma::vec time_unique(y.n_rows);

  // use sorted y times to count the number of unique.
  // also define number at risk as the sum of the weights
  arma::uword n_unique = 1; // set at 1 to account for the first time
  n_risk = weights[0];      // see above
  time_unique(0) = y(0,0);  // see above

  for(person = 1; person < y.n_rows; person++){

    if(y(person-1, 0) != y(person,0)){
      time_unique(n_unique) = y(person,0);
      n_unique++;
    }

    n_risk += weights[person];

  }

  // drop the extra zeros from time_unique
  time_unique = time_unique(arma::span(0, n_unique-1));
  //Rcout << n_risk << std::endl;
  //Rcout << n_unique << std::endl;

  // reset for next loop
  person = 0;
  km_counter = 0;
  double km = 1.0;
  double person_time;
  arma::vec kmvec(n_unique);

  do{

    person_time = y(person, 0);
    n_dead = 0;
    n_cens = 0;
    n_risk_sub = 0;

    while(y(person, 0) == person_time){

      n_risk_sub += weights[person];

      if(y(person, 1) == 1){
        n_dead += weights[person];
      } else {
        n_cens += weights[person];
      }

      if(person == y.n_rows-1) break;
      person++;

    }

    //Rcout << "n_risk: " << n_risk << std::endl;
    //Rcout << "n_dead: " << n_dead << std::endl;
    //Rcout << "n_risk: " << n_risk << std::endl;

    km = km * (n_risk - n_dead) / n_risk;

    //Rcout << "km: " << km << std::endl;

    kmvec[km_counter] = km;

    n_risk -= n_risk_sub;

    km_counter++;

  } while (km_counter < n_unique);

  //Rcout << time_unique.n_rows << std::endl;
  //Rcout << kmvec.n_rows << std::endl;

  return(arma::join_horiz(time_unique, kmvec));

}

// [[Rcpp::export]]
arma::mat node_summarize(arma::mat& y,
                         arma::uvec& node_assignments,
                         arma::uvec& weights,
                         arma::uword& nodes_max){


  // allocate memory for output
  arma::mat out(nodes_max, 2);
  arma::uword row_index;

  // loop through the matrix, once, by row
  for(arma::uword i = 0; i < node_assignments.size(); i++){

    // subtract 1 from current value of nodes to align
    // with c++ index starting at 0.

    // add the current event value from i'th row of Y
    // to the current bucket in the output, which
    // is determined by the current value of nodes

    row_index = node_assignments[i]-1;

    out(row_index, 0) += y(i, 1) * weights(i);
    out(row_index, 1) += weights(i);

  }

  return(out);

}


// [[Rcpp::export]]
void find_cutpoints(arma::vec& cp,
                    const arma::vec& lc,
                    const arma::mat& y,
                    const arma::uword& leaf_min_obs,
                    const arma::uword& leaf_min_events){

  arma::vec lc_uni;
  arma::uword lc_size = lc.size();
  arma::uword cp_size = cp.size();
  arma::uword i, j;
  arma::uvec index;
  arma::uvec lc_index;

  double n_events_left, n_events_right, n_obs_left, n_obs_right;

  for(i = 0; i < cp_size; i++){
    cp[i] = R_PosInf;
  }

  if(lc_size > 20){

    // get unique values in the head and tail, using the rough sort:
    // - greater values tend to be at the head of x
    // - smaller values tend to be at the tail of x
    // This is a shallow way to check the number of unique values
    lc_uni = arma::unique(
      arma::join_cols(lc.head(10),
                      lc.tail(10))
    );


  } else {

    lc_uni = arma::unique(lc);

  }

  // note that arma::unique() sorts the output from low to high

  if(lc_uni.size() == 1){

    // if there is only one unique value,
    // then there is no cut-point to be made.
    // --> do nothing

  } else if (lc_uni.size() <= cp_size){
    // if there are k <= n_cutpoint unique values:
    // this means we have k-1 valid cutpoints
    // example: lc has unique values of -1, 0, and 1.
    // if -1 is the cut-point, then all values <= -1 go left, >-1 go right
    // if 0 is the cut-point, then all values <= 0 go left, > 0 go right
    // the max (1 in this case) is not a valid cutpoint, bc nothing is > max
    // each unique value of lc is a valid cut-point.

    //Rcout << "lc_uni: " << lc_uni.t() << std::endl;

    for(i = 0; i < lc_uni.size(); i++){

      n_events_left = 0;
      n_events_right = 0;
      n_obs_left = 0;
      n_obs_right = 0;

      //Rcout << "lc: " << lc[i] << std::endl;

      for(j = 0; j < lc_size; j++){

        if(lc[j] <= lc_uni[i]){
          n_events_left += y(j, 1);
          n_obs_left++;
        } else {
          n_events_right += y(j, 1);
          n_obs_right++;
        }

        // Rcout << "n_events_left: " << n_events_left << std::endl;
        // Rcout << "n_events_right: " << n_events_right << std::endl;
        // Rcout << "n_obs_left: " << n_obs_left << std::endl;
        // Rcout << "n_obs_right: " << n_obs_right << std::endl;

        // check if left node has enough events and observations
        // (don't forget to make sure the right node has enough as well)
        if(n_events_left >= leaf_min_events &&
           n_events_right >= leaf_min_events &&
           n_obs_left >= leaf_min_obs &&
           n_obs_right >= leaf_min_obs){
          cp[i] = lc_uni[i];
          break;
        }

      }

    }

  } else { // if there are k > n_cutpoint unique values
    // assume this is a continuous linear combination vector

    arma::vec probs;

    if(cp_size >= 3){
      probs = arma::linspace<arma::vec>(0.20, 0.80, cp_size);
    } else if (cp_size == 2){
      probs = {0.333, 0.666};
    } else {
      probs = {0.50};
    }

    arma::vec lc_quants;

    if(lc_size >= 100){
      lc_quants = arma::quantile(
        lc(arma::linspace<arma::uvec>(0, lc_size-1, 100)),
        probs
      );
    } else {
      lc_quants = arma::quantile(lc, probs);
    }


    // arma::vec lc_quants_true = arma::quantile(lc, probs);
    //
    // Rcout << "approx q: " << lc_quants.t() << std::endl;
    //
    // Rcout << "true q: " << lc_quants_true.t() << std::endl;


    //Rcout << "lc: " << lc[i] << std::endl;
    for(i = 0; i < cp_size; i++){

      n_events_left = 0;
      n_events_right = 0;
      n_obs_left = 0;
      n_obs_right = 0;

      for(j = 0; j < lc_size; j++){

        if(lc[j] <= lc_quants[i]){
          n_events_left += y(j, 1);
          n_obs_left++;
        } else {
          n_events_right += y(j, 1);
          n_obs_right++;
        }

        // check if left node has enough events and observations
        if(n_events_left >= leaf_min_events &&
           n_events_right >= leaf_min_events &&
           n_obs_left >= leaf_min_obs &&
           n_obs_right >= leaf_min_obs){
          cp[i] = lc_quants[i];

          // Rcout << "n_events_left: " << n_events_left << std::endl;
          // Rcout << "n_events_right: " << n_events_right << std::endl;
          // Rcout << "n_obs_left: " << n_obs_left << std::endl;
          // Rcout << "n_obs_right: " << n_obs_right << std::endl;

          break;
        }

      }

    }

  }

  //Rcout << cp.t() << std::endl;

}

// [[Rcpp::export]]
double log_rank_test(arma::mat& y,
                     arma::vec& g){

  arma::uword n = y.n_rows;
  double Y = n;
  double Y1 = arma::sum(g);

  arma::uword lwr = 0;
  arma::uword upr = 0;
  arma::uword count = 1; // starts at 1 to mimic size
  arma::uword i;

  for(i=0; i<n; i++){

    if(y(i, 1) == 0){

      upr++;
      count++;

    } else {

      break;

    }

  }

  double d = 0;
  double d1 = 0;

  for(i = lwr; i <= upr; i++){
    d += y(i, 1);
    d1 += y(i, 1) * g(i);
  }

  double e1 = Y1 * d / Y;
  double e0 = (Y - Y1) * d / Y;
  double o1 = d1;
  double o0 = d - d1;

  double V = (Y-Y1) * Y1 * d * (Y-d) / (pow(Y,2) * (Y-1));

  Y -= count;

  for(i = lwr; i <= upr; i++){
    Y1 -= g(i);
  }

  lwr=upr+1;

  for( ; ; ){

    // Rcout << "e1: " << e1 << "; ";
    // Rcout << "e0: " << e0 << "; ";
    // Rcout << "o1: " << o1 << "; ";
    // Rcout << "o0: " << o0 << "; ";
    // Rcout << "Y1: " << Y1 << "; ";
    // Rcout << "Y: " << Y << "; ";
    // Rcout << std::endl;

    upr = lwr;
    count = 1;


    while( (y(upr, 1) == 0) & (upr < n-1) ){
      upr++;
      count++;
    }

    if(upr==n-1){
      if( y(upr, 1) == 0 ){
        break;
      } else {

        d = 0; d1 = 0;

        for(i = lwr; i <= upr; i++){
          d += y(i, 1);
          d1 += y(i, 1) * g(i);
        }

        e1 += (Y1 * d/Y);
        e0 += ((Y-Y1) * d/Y);
        o1 += d1;
        o0 += d-d1;

        V += (Y-Y1) * Y1 * d * (Y-d) / (pow(Y, 2) * (Y-1));
        Y -= count;

        for(i = lwr; i <= upr; i++) Y1 -= g(i);

        break;
      }
    }

    d = 0; d1 = 0;

    for(i = lwr; i <= upr; i++){
      d += y(i, 1);
      d1 += y(i, 1) * g(i);
    }

    e1 += (Y1*d / Y);
    e0 += ((Y-Y1) * d/Y);
    o1 += d1;
    o0 += d - d1;

    V += (Y-Y1) * Y1 * d * (Y-d) / (pow(Y,2) * (Y-1));
    Y -= count;

    for(i = lwr; i <= upr; i++) Y1 -= g(i);

    lwr=upr+1;

    if(Y==1) break;

  }

  return pow(o1-e1,2) / V;

}

// TODO: write x_scale function (copy from rcpp version)
// [[Rcpp::export]]
NumericMatrix x_scale(NumericMatrix x_mat,
                      IntegerVector x_wts,
                      LogicalVector do_scale){

  int ncol = x_mat.ncol();
  int nrow = x_mat.nrow();
  int i, person;
  double x_scale_i, x_mean_i;
  // set aside memory for outputs
  // first column holds the mean values
  // second column holds the scale values
  NumericMatrix out(ncol, 2);
  NumericMatrix::Column means = out( _ , 0);   // Reference to column 1
  NumericMatrix::Column scales = out( _ , 1);  // Reference to column 2

  int x_wts_sum = sum(x_wts);

  for(i = 0; i < ncol; i++) {

    if (do_scale[i] == false) {

      scales[i] = 1.0;
      means[i] = 0.0;

    } else {

      x_mean_i = 0.0;

      for (person = 0; person < nrow; person++){
        x_mean_i = x_mean_i + x_wts[person] * x_mat(person, i);
      }

      means[i] = x_mean_i / x_wts_sum;

      for (person = 0; person < nrow; person++){
        x_mat(person, i) -= means[i];
      }

      x_scale_i = 0.0;

      for (person = 0; person < nrow; person++){
        x_scale_i = x_scale_i + x_wts[person] * abs(x_mat(person, i));
      }

      if(x_scale_i > 0)
        scales[i] = x_wts_sum / x_scale_i;
      else
        scales[i] = 1.0; // rare case of constant covariate;

      for (person = 0; person < nrow; person++){
        x_mat(person, i) *= scales[i];
      }

    }

  }

  return(out);

}

// [[Rcpp::export]]
void cholesky(arma::mat& matrix){

  double eps = 0;
  double toler = 1e-8;
  double pivot, temp;

  arma::uword i, j, k;
  arma::uword n = matrix.n_cols;

  for(i = 0; i < n; i++){

    if(matrix(i,i) > eps) eps = matrix(i,i);

    // copy upper right values to bottom left
    for(j = (i+1); j<n; j++){
      matrix(j,i) = matrix(i,j);
    }
  }

  if (eps == 0)
    eps = toler; // no positive diagonals!
  else
    eps = eps * toler;

  for (i = 0; i < n; i++) {

    pivot = matrix(i, i);

    if (pivot == std::numeric_limits<double>::infinity() || pivot < eps) {

      matrix(i, i) = 0;

    } else {

      for(j = (i+1); j < n; j++){

        temp = matrix(j,i) / pivot;
        matrix(j,i) = temp;
        matrix(j,j) -= temp*temp*pivot;

        for(k = (j+1); k < n; k++){

          matrix(k, j) -= temp * matrix(k, i);

        }

      }

    }

  }

}

// [[Rcpp::export]]
void cholesky_solve(arma::mat& matrix,
                    arma::vec& y){

  arma::uword n = matrix.n_cols;
  arma::uword i, j;
  double temp;

  for (i = 0; i < n; i++) {

    // Rcout << i << std::endl;

    temp = y[i];

    for (j = 0; j < i; j++){

      // Rcout << " " << j << std::endl;

      temp -= y[j] * matrix(i, j);
      y[i] = temp;

    }

  }

  // A hack so that i is never < 0 (out of index for uwords)
  arma::uword ii;

  for (i = (n); i >= 1; i--){

    // Rcout << i << std::endl;
    ii = i-1;

    if (matrix(ii, ii) == 0){

      y[ii] =0;

    } else {

      temp = y[ii] / matrix(ii, ii);

      for (j = ii+1; j < n; j++){
        // Rcout << " " << j << std::endl;
        temp -= y[j] * matrix(j, ii);
      }

      y[ii] = temp;

    }

  }

}

// [[Rcpp::export]]
void cholesky_invert(arma::mat& matrix){

  double temp;
  arma::uword i,j,k;
  arma::uword n = matrix.n_cols;

  /*
   ** invert the cholesky in the lower triangle
   **   take full advantage of the cholesky's diagonal of 1's
   */
  for (i=0; i<n; i++){

    if (matrix(i,i) >0) {

      matrix(i,i) = 1.0 / matrix(i,i);

      for (j=(i+1); j<n; j++) {

        matrix(j, i) = -matrix(j, i);

        for (k=0; k<i; k++){
          matrix(j, k) += matrix(j, i) * matrix(i, k);
        }

      }

    }

  }

  /*
   ** lower triangle now contains inverse of cholesky
   ** calculate F'DF (inverse of cholesky decomp process) to get inverse
   **   of original matrix
   */
  for (i=0; i<n; i++) {

    if (matrix(i, i) == 0) {

      for (j=0; j<i; j++) matrix(i, j) = 0;
      for (j=i; j<n; j++) matrix(j, i) = 0;

    } else {

      for (j=(i+1); j<n; j++) {

        temp = matrix(j, i) * matrix(j, j);

        if (j!=i) matrix(i, j) = temp;

        for (k=i; k<j; k++){
          matrix(i, k) += temp*matrix(j, k);
        }

      }

    }

  }

}


// [[Rcpp::export]]
void newtraph_cph_one_iter (const arma::mat& x,
                            const arma::mat& y,
                            const arma::uvec& weights,
                            arma::vec& u,
                            arma::vec& a,
                            arma::vec& a2,
                            arma::mat& imat,
                            arma::mat& cmat,
                            arma::mat& cmat2,
                            int method){


  double  wtave;
  double  temp1, temp2;
  double  person_time;

  arma::uword deadwt;

  arma::uword i, j, k, ndead, denom2;

  arma::uword nrisk = 0;
  arma::uword denom = 0;

  arma::uword person = x.n_rows - 1;
  arma::uword nvar = x.n_cols;

  bool break_loop = false;
  // this loop has a strange break condition to accomodate
  // the restriction that a uvec (uword) cannot be < 0
  for ( ; ; ){

    // Rcout << "- person: "    << person;
    // Rcout << "; y_status: "  << y(person,1);
    // Rcout << "; y_time: "    << y(person,0);
    // Rcout << "; u: "         << u.t();
    // Rcout << std::endl;

    person_time = y(person, 0); // time of event for current person
    ndead  = 0 ; // number of deaths at this time point
    deadwt = 0 ; // sum of weights for the deaths
    denom2 = 0 ; // sum of weighted risks for the deaths

    // walk through this set of tied times
    while(y(person, 0) == person_time){

      nrisk++;

      if (y(person, 1) == 0) {

        denom += weights[person];

        /* a contains weighted sums of x, cmat sums of squares */

        for (i=0; i<nvar; i++) {

          temp1 = weights[person] * x(person, i);

          a[i] += temp1;

          for (j=0; j<=i; j++){
            cmat(j, i) += temp1 * x(person, j);
          }


        }

      } else {

        ndead++;

        deadwt += weights[person];
        denom2 += weights[person];

        for (i=0; i<nvar; i++) {

          temp1 = weights[person] * x(person, i);

          u[i]  += temp1;
          a2[i] += temp1;

          for (j=0; j<=i; j++){
            cmat2(j, i) += temp1 * x(person, j);
          }

        }

      }

      if(person == 0){
        break_loop = true;
        break;
      }
      person--;

    }

    //Rcout << "; ndead: "         << ndead;
    //Rcout << std::endl;

    //Rcout << "imat: " << std::endl << imat << std::endl;

    // we need to add to the main terms
    if (ndead > 0) {

      if (method == 0 || ndead == 1) { // Breslow

        denom  += denom2;

        for (i=0; i<nvar; i++) {

          a[i]  += a2[i];
          temp2  = a[i] / denom;
          u[i]  -=  deadwt * temp2;

          for (j=0; j<=i; j++) {
            cmat(j, i) += cmat2(j, i);
            imat(j, i) += deadwt * (cmat(j, i) - temp2 * a[j]) / denom;
          }

        }

      } else {
        /* Efron
         **  If there are 3 deaths we have 3 terms: in the first the
         **  three deaths are all in, in the second they are 2/3
         **  in the sums, and in the last 1/3 in the sum.  Let k go
         **  1 to ndead: we sequentially add a2/ndead and cmat2/ndead
         **  and efron_wt/ndead to the totals.
         */
        wtave = deadwt/ndead;

        for (k=0; k<ndead; k++) {

          denom  += denom2 / ndead;

          for (i=0; i<nvar; i++) {

            a[i] += a2[i] / ndead;
            temp2 = a[i]  / denom;
            u[i] -= wtave * temp2;

            for (j=0; j<=i; j++) {
              cmat(j, i) += cmat2(j, i) / ndead;
              imat(j, i) += wtave * (cmat(j, i) - temp2 * a[j]) / denom;
            }

          }

        }

      }

      for (i=0; i<nvar; i++) {

        a2[i]=0;
        for (j=0; j<nvar; j++) cmat2(j,i)=0;

      }

    }

    if(break_loop == true) break;
    if(person == 0) break_loop = true;

  }

  //Rcout << "- u final: " << u.t() << std::endl;

}


// [[Rcpp::export]]
String make_node_name(const arma::uword& part){

  std::ostringstream oss;
  oss << "node_" << part;
  String out = oss.str();

  return out;
}

// [[Rcpp::export]]
bool any_cps_valid(arma::vec& x){

  for(arma::uword i = 0; i < x.size(); i++){
    if (x[i] < R_PosInf) return(true);
  }

  return false;

}

// [[Rcpp::export]]
List ostree_fit_arma(arma::mat& x,
                     arma::mat& y,
                     const arma::uword& mtry = 4,
                     const arma::uword& n_cps = 5,
                     const arma::uword& leaf_min_events = 5,
                     const arma::uword& leaf_min_obs = 10,
                     const bool& verbose = false){

  int n_obs = x.n_rows;
  int n_col = x.n_cols;
  double prob_sampled = 1.0/n_obs;

  if(verbose == true){
    Rcout << std::endl << "N obs total: " << n_obs << std::endl;
    Rcout << "N columns total: " << n_col << std::endl << std::endl;
  }

  // ----------------------------------------------------
  // ---- sample weights to mimic a bootstrap sample ----
  // ----------------------------------------------------

  // s is the number of times you might get selected into
  // a bootstrap sample. Realistically this won't be >10.
  IntegerVector s = seq(0,10);

  // compute probability of being selected into the bootstrap
  // 0 times, 1, times, ..., 9 times, or 10 times.
  NumericVector probs = dbinom(s, n_obs, prob_sampled, false);

  // sampling with replacement
  arma::uvec weights = as<arma::uvec>(sample(s, n_obs, true, probs));

  if(verbose == true){

    Rcout << "weight vector (first 7): " << std::endl <<
      weights.head(7).t() << std::endl;

    Rcout << "weight vector (last 7): " << std::endl <<
      weights.tail(7).t() << std::endl;

  }

  // vector indicating which rows are in the bootstrap sample
  // note: obs i is in the sample if weights[i] > 0
  arma::uvec rows_inbag = arma::find(weights);

  // drop the zero's from weights (weights only contain in-bag data now)
  weights = weights(rows_inbag);
  arma::uword n_obs_wts = sum(weights);

  // sub-matrices containing just the bootstrapped rows
  // used instead of the original y and x matrix
  // (Question: is this better than using just a row index
  //            and keeping the original x and y matrices?)
  arma::mat y_inbag = y.rows(rows_inbag);
  arma::mat x_inbag = x.rows(rows_inbag);


  if(verbose == true){

    Rcout << "N obs in bootstrap aggregate (in-bag) sample: " <<
      rows_inbag.size() << std::endl;

    Rcout << "max times one obs sampled: " <<
      arma::max(weights) << std::endl;

    Rcout << "sum of bootstrap weights: " << n_obs_wts << std::endl;

  }

  // ----------------------------------------
  // ---- preallocate memory for outputs ----
  // ----------------------------------------

  // the maximum number of nodes possible given the tree params
  // TODO: Underestimate max nodes and add memory if needed.
  arma::uword max_nodes = 2 * std::ceil(n_obs_wts / leaf_min_events) - 1;
  // node counter tracks progression through the allocated memory
  // arma::uword node_counter = 0;
  // beta coefficients for linear combinations
  arma::mat betas(mtry, max_nodes);
  // column indices to pair with beta coefficients
  arma::umat col_indices(mtry, max_nodes);
  // cutpoints for linear combination splitting
  arma::vec cutpoints(max_nodes);
  // child nodes to progress through the tree
  arma::uvec children_left(max_nodes),  children_right(max_nodes);

  // --------------------------------------------
  // ---- initialize parameters to grow tree ----
  // --------------------------------------------

  // the nodes vector keeps track of which node each observation is in.
  // everyone starts at node 0, the root of the tree;
  arma::uvec node_assignments( rows_inbag.size() );

  // nodes_max indicates the current highest value of nodes;
  // this determines what new nodes will be labeled as.
  arma::uword nodes_max = 0;

  // nodes_to_grow
  // - indicates nodes with enough obs and events to be split
  // - updates after every iteration
  // - iterations end if this vector is empty.
  arma::uvec nodes_to_grow {0};

  // sub-matrices containing just the rows in the current node
  // used to count events in nodes and make probs in leaf nodes
  arma::mat y_node, x_node;
  // weights will also be subsetted here.
  arma::uvec weights_node;
  // vectors to make x_node, y_node, and weights_node
  arma::uvec rows_node, cols_node;

  // iterators for the main loop
  arma::uword i, j;
  arma::uvec::iterator node, person;

  // dealing with sampling of columns
  arma::uword n_cols_to_sample;
  // columns to be sampled
  arma::uvec cols_to_sample_01(n_col);
  // if the number of cols to sample is < mtry,
  // we will set mtry_temp = n_cols_to_sample
  arma::uword mtry_temp;
  // columns to sample
  arma::uvec cols_to_sample;

  // vectors and matrices for newton raphson
  arma::vec u(mtry);
  arma::vec a(mtry);
  arma::vec a2(mtry);
  arma::mat imat(mtry, mtry);
  arma::mat cmat(mtry, mtry);
  arma::mat cmat2(mtry, mtry);

  // vectors from solving newton raphson (first iter)
  arma::vec beta, lc, lc_temp;

  // vector to hold cut-points
  arma::vec cp(n_cps);

  // data to check event/sample size in each part
  arma::mat node_summary;

  // container for leaf node data
  List leaf_nodes;
  //each leaf node is a matrix with varying row count
  arma::mat leaf;

  // subset data for the current node
  arma::mat y_leaf;
  arma::uvec rows_leaf;
  arma::uvec weights_leaf;
  arma::uword nn_left;
  arma::uword nn_right;
  String node_name;

  // data from nodes that were grown in current iteration
  arma::mat y_grown;
  arma::uvec nodes_grown;
  arma::uvec weights_grown;

  // ----------------------
  // ---- main do loop ----
  // ----------------------

  do { // looping through nodes that need to be split

    arma::uvec rows_node_combined;

    for(node = nodes_to_grow.begin(); node != nodes_to_grow.end(); ++node){

      if(verbose == true){
        Rcout << std::endl <<
          "---- Growing node " << *node << " ----" <<
            std::endl << std::endl;
      }

      mtry_temp = mtry;

      u.fill(0);
      a.fill(0);
      a2.fill(0);
      imat.fill(0);
      cmat.fill(0);
      cmat2.fill(0);
      cols_to_sample_01.fill(0);

      if(nodes_to_grow[0] == 0){

        rows_node = arma::linspace<arma::uvec>(0,
                                               rows_inbag.size()-1,
                                               rows_inbag.size());

      } else {

        rows_node = arma::find(node_assignments == *node);

      }

      // check for constant columns
      // constant means constant within the rows where events occurred

      if(verbose == true){
        Rcout <<
          "empty cols_to_sample: " <<
            cols_to_sample_01.t() <<
              std::endl;
      }



      for(j = 0; j < cols_to_sample_01.size(); j++){

        double reference = R_PosInf;

        for(person = rows_node.begin()+1; person != rows_node.end(); ++person){

          if(y_inbag(*person, 1) == 1){

            if (reference < R_PosInf){

              if(x_inbag(*person, j) != reference){

                cols_to_sample_01[j] = 1;
                //Rcout << "column " << k << " is okay" << std::endl;
                break;

              }

            } else {

              reference = x_inbag(*person, j);

            }


          }

        }

      }

      if(verbose == true){
        Rcout <<
          "full  cols_to_sample: " <<
            cols_to_sample_01.t() <<
              std::endl;
      }

      n_cols_to_sample = arma::sum(cols_to_sample_01);

      if(n_cols_to_sample < mtry) mtry_temp = n_cols_to_sample;

      if(verbose == true && mtry_temp < mtry){
        Rcout << "mtry_temp: " << mtry_temp << std::endl;
      }


      // a hack to convert mtry_temp to an int
      int mtry_temp_int = 0;

      for(i = 0; i < mtry_temp; i++) mtry_temp_int++;

      //Rcout << "mtry_temp_int: " << mtry_temp_int << std::endl;

      cols_to_sample = arma::find(cols_to_sample_01);

      cols_node = Rcpp::RcppArmadillo::sample(cols_to_sample,
                                              mtry_temp_int,
                                              false);

      if(verbose == true){
        Rcout << "cols sampled: " << cols_node.t() << std::endl;
      }


      x_node = x_inbag(rows_node, cols_node);
      y_node = y_inbag.rows(rows_node);
      weights_node = weights(rows_node);

      // Rcout << "x_node: " << std::endl << x_node << std::endl;
      //
      // Rcout << "y_node: " << std::endl << y_node << std::endl;

      //Rcout << "weights: " << weights.t() << std::endl;

      // Rcout << "u: " << u.t() << std::endl;

      newtraph_cph_one_iter(x_node, y_node, weights_node,
                            u, a, a2, imat, cmat, cmat2, 0);

      //Rcout << "imat: " << std::endl << imat << std::endl;
      cholesky(imat);
      //Rcout << "imat after cholesky: " << std::endl << imat << std::endl;
      cholesky_solve(imat, u);

      arma::vec beta(mtry_temp);

      for(i = 0; i < mtry_temp; i++) beta[i] = u[i];

      if(verbose == true){
        Rcout << "beta: " << beta.t() << std::endl;
      }

      lc = x_node * beta;

      find_cutpoints(cp, lc, y_node, leaf_min_obs, leaf_min_events);

      // if there is at least one valid cut-point, we can grow
      // two more nodes from the current node.

      if(any_cps_valid(cp) == true){

        rows_node_combined = arma::join_cols(rows_node_combined,
                                             rows_node);

        if(verbose == true){
          Rcout << "cut points: " << cp.t() << std::endl;
        }

        arma::vec g(lc.size());

        double lrstat_max = 0;
        double lrstat, cp_max;

        for(i = 0; i < cp.size(); i++){

          g.fill(0);

          if(cp[i] < R_PosInf){

            for(j = 0; j < lc.size(); j++){
              if(lc[j] > cp[i]) g[j] = 1;
            }

            if(verbose == true){
              Rcout << "number of g == 1: " <<
                arma::sum(g) <<
                  std::endl;
              Rcout <<
                "number of g == 0: " <<
                  g.size() - arma::sum(g) <<
                    std::endl;
            }

            lrstat = log_rank_test(y_node, g);

            if(verbose == true){
              Rcout << lrstat << " with cut-point " << cp[i];
            }


            if(lrstat > lrstat_max){
              lrstat_max = lrstat;
              cp_max = cp[i];
              if(verbose == true) Rcout << ", a new max!";
            }

            if(verbose == true) Rcout << std::endl;

          }

        }

        nn_left   = nodes_max + 1;
        nn_right  = nodes_max + 2;
        nodes_max = nn_right;

        i = 0;

        for(person = rows_node.begin(); person != rows_node.end(); ++person){

          if(lc[i] <= cp_max){

            node_assignments[*person] = nn_left;

          } else {

            node_assignments[*person] = nn_right;

          }

          i++;

        }

        for(i = 0; i < mtry_temp; i++){
          betas(i, *node) = beta(i);
          col_indices(i, *node) = cols_node(i);
        }

        children_left[*node] = nn_left;
        children_right[*node] = nn_right;
        cutpoints[*node] = cp_max;

        //node_counter++;


      } else {

        // some nodes can make the loop go a longer time because they
        // technically have enough obs and events to be split but
        // its hard to find a cut-point that can meet minimum number
        // of obs and events in the two nodes that would be created
        // from the split. So we blacklist nodes that do that.

        if(verbose == true){
          Rcout << "dont grow: " << *node << std::endl;
        }

        leaf = leaf_surv_small(y_node, weights_node);

        if(verbose == true){
          Rcout << "creating a new leaf: node_" << *node << std::endl;
        }

        leaf_nodes[make_node_name(*node)] = leaf;

        //node_counter++;

      }

    }

    y_grown       = y_inbag.rows(rows_node_combined);
    nodes_grown   = node_assignments(rows_node_combined);
    weights_grown = weights(rows_node_combined);

    node_summary = node_summarize(y_grown,
                                  nodes_grown,
                                  weights_grown,
                                  nodes_max);

    if(verbose == true){
      Rcout << std::endl;
      Rcout << "events by part: " << std::endl;
      Rcout << node_summary << std::endl;
    }

    arma::uvec nodes_to_grow_temp(node_summary.n_rows);

    for(i = 0; i < node_summary.n_rows; i++){

      if(node_summary(i, 0) >= 2 * leaf_min_events &&
         node_summary(i, 1) >= 2 * leaf_min_obs){
        // split it

        // use i+1; the ith row of node_summary is the i+1 node
        nodes_to_grow_temp(i) = i + 1;

      } else if (node_summary(i, 0) > 0 &&
                 node_summary(i, 1) > 0) {

        // a new leaf
        // use i+1; nodes starts at 1 and i starts at 0
        rows_leaf    = arma::find(node_assignments == i+1);
        y_leaf       = y_inbag.rows(rows_leaf);
        weights_leaf = weights(rows_leaf);
        leaf     = leaf_surv_small(y_leaf, weights_leaf);

        if(verbose == true){
          Rcout << "created new leaf: node_" << i+1 << std::endl;
        }

        leaf_nodes[make_node_name(i+1)] = leaf;
        //node_counter++;

      }

    }

    nodes_to_grow = nodes_to_grow_temp(arma::find(nodes_to_grow_temp));

    if(verbose == true){
      Rcout << "nodes to grow: " << nodes_to_grow.t() << std::endl;
    }

  } while (nodes_to_grow.size() > 0);

  return(
    List::create(
      Named("leaves") = leaf_nodes,
      _["betas"] = betas.cols(arma::span(0, nodes_max)),
      _["col_indices"] = col_indices.cols(arma::span(0, nodes_max)),
      _["cut_points"] = cutpoints(arma::span(0, nodes_max)),
      _["children_left"] = children_left(arma::span(0, nodes_max)),
      _["children_right"] = children_right(arma::span(0, nodes_max)),
      _["mtry"] = mtry
    )
  );




}


// [[Rcpp::export]]
arma::uvec ostree_pred_leaf1(const arma::mat& x_new,
                             const arma::mat& betas,
                             const arma::umat& col_indices,
                             const arma::vec& cut_points,
                             const arma::vec& children_left,
                             const arma::vec& children_right){

  // allocate memory for output
  arma::uvec out(x_new.n_rows);
  arma::uword i, j, k;
  arma::uword n_nodes = betas.n_cols;
  arma::uvec obs_in_node;

  double temp;

  for(i = 0; i < n_nodes; i++){

    //Rcout << "i: " << i << std::endl;

    if(children_left(i) != 0){

      obs_in_node = arma::find(out == i);

      if(obs_in_node.size() > 0){

        for(j = 0; j < obs_in_node.size(); j++){

          //Rcout << " - j: " << j << std::endl;
          temp = 0;

          for(k = 0; k < betas.n_rows; k++){

            //Rcout << " -- k: " << k << std::endl;

            temp += x_new( obs_in_node(j), col_indices(k, i) ) * betas(k, i);

          }


          if(temp <= cut_points[i]) {

            out(obs_in_node(j)) = children_left(i);

          } else {

            out(obs_in_node(j)) = children_right(i);

          }

        }

      }

    }

  }

  return(out);

}

// [[Rcpp::export]]
arma::uvec ostree_pred_leaf2(const arma::mat& x_new,
                             const arma::mat& betas,
                             const arma::umat& col_indices,
                             const arma::vec& cut_points,
                             const arma::vec& children_left,
                             const arma::vec& children_right){

  // allocate memory for output
  arma::uvec out(x_new.n_rows);
  arma::uword i, k;
  arma::uvec::iterator j;
  arma::uword n_nodes = betas.n_cols;
  arma::uvec obs_in_node;

  double temp;

  for(i = 0; i < n_nodes; i++){

    //Rcout << "i: " << i << std::endl;

    if(children_left(i) != 0){

      obs_in_node = arma::find(out == i);

      if(obs_in_node.size() > 0){

        for(j = obs_in_node.begin(); j != obs_in_node.end(); ++j){

          //Rcout << " - j: " << j << std::endl;
          temp = 0;

          for(k = 0; k < betas.n_rows; k++){

            //Rcout << " -- k: " << k << std::endl;

            temp += x_new( *j, col_indices(k, i) ) * betas(k, i);

          }


          if(temp <= cut_points[i]) {

            out(*j) = children_left(i);

          } else {

            out(*j) = children_right(i);

          }

        }

      }

    }

  }

  return(out);

}

// [[Rcpp::export]]
arma::uvec ostree_pred_leaf3(const arma::mat& x_new,
                             const arma::mat& betas,
                             const arma::umat& col_indices,
                             const arma::vec& cut_points,
                             const arma::vec& children_left,
                             const arma::vec& children_right){

  // allocate memory for output
  arma::uvec out(x_new.n_rows);
  arma::uword i, k;
  arma::uword n_nodes = betas.n_cols;
  arma::uvec obs_in_node;

  arma::mat x_node;
  arma::vec lc;

  for(i = 0; i < n_nodes; i++){

    //Rcout << "i: " << i << std::endl;

    if(children_left(i) != 0){

      obs_in_node = arma::find(out == i);

      if(obs_in_node.size() > 0){

        x_node = x_new(obs_in_node, col_indices.col(i));

        lc = x_node * betas.col(i);

        for(k = 0; k < lc.size(); k++)

          if(lc(k) <= cut_points(i)) {

            out(obs_in_node(k)) = children_left(i);

          } else {

            out(obs_in_node(k)) = children_right(i);

          }

      }

    }

  }

  return(out);

}

// [[Rcpp::export]]
arma::uvec ostree_pred_leaf4(const arma::mat& x_new,
                             const arma::mat& betas,
                             const arma::umat& col_indices,
                             const arma::vec& cut_points,
                             const arma::vec& children_left,
                             const arma::vec& children_right){

  // allocate memory for output
  arma::uvec out(x_new.n_rows);
  arma::uword i, j, k;
  arma::uword n_nodes = betas.n_cols;
  arma::uvec obs_in_node;

  arma::vec lc;

  for(i = 0; i < n_nodes; i++){

    //Rcout << "i: " << i << std::endl;

    if(children_left(i) != 0){

      obs_in_node = arma::find(out == i);

      if(obs_in_node.size() > 0){

        lc = x_new(obs_in_node, col_indices.col(i)) * betas.col(i);

        for(j = 0; j < lc.size(); j++){

          k = obs_in_node(j);

          if(lc(j) <= cut_points(i)) {

            out(k) = children_left(i);

          } else {

            out(k) = children_right(i);

          }

        }

      }

    }

  }

  return(out);

}




