#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat x_mean_sd(const arma::mat& x){

  // 0 indicates statistic is computed by column
  arma::mat x_mean = arma::mean(x, 0);
  arma::mat x_sd = arma::stddev(x, 0, 0);
  return(arma::join_vert(x_mean, x_sd));

}

// [[Rcpp::export]]
void x_scale(arma::mat& x,
             arma::mat& mean_sd){

  for(arma::uword i = 0; i < x.n_cols; i++){
    x.unsafe_col(i) = (x.unsafe_col(i) - mean_sd(0,i)) / mean_sd(1,i);
  }

}


// [[Rcpp::export]]
arma::mat x_scale_cph(arma::mat& x_mat,
                      arma::uvec& weights){

  // set aside memory for outputs
  // first column holds the mean values
  // second column holds the scale values

  arma::mat out(x_mat.n_cols, 2);
  arma::vec means = out.unsafe_col(0);   // Reference to column 1
  arma::vec scales = out.unsafe_col(1);  // Reference to column 2
  arma::uword weights_sum = arma::sum(weights);

  for(arma::uword i = 0; i < x_mat.n_cols; i++) {

    arma::vec x_i = x_mat.unsafe_col(i);

    means(i) = arma::sum( weights % x_i ) / weights_sum;

    x_i -= means(i);

    scales(i) = arma::sum(weights % arma::abs(x_i));

    if(scales(i) > 0)
      scales(i) = weights_sum / scales(i);
    else
      scales(i) = 1.0; // rare case of constant covariate;

    x_i *= scales(i);

  }

  return(out);

}

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
void find_cutpoints_ctns(arma::vec& cp,
                         const arma::vec& lc,
                         const arma::mat& y,
                         const arma::uvec weights,
                         const arma::uword& leaf_min_obs,
                         const arma::uword& leaf_min_events){

  // sort the lc vector
  arma::uvec lc_sort = arma::sort_index(lc);
  // iterate over y and weights; low to high lc and vice versa.
  arma::uvec::iterator it;
  // track events and obs to determine cutpoint validity.
  // keep track of minimal and maximal valid cut-points.
  double n_events, n_obs, cp_min, cp_max;

  // beginning from lowest values,
  // cycle up until you have leaf_min_obs and leaf_min_events
  // designate the value of lc at this spot as the minimal cutpoint

  for (it = lc_sort.begin(); it != lc_sort.end(); ++it){

    n_events += y(*it, 1) * weights(*it);
    n_obs += weights(*it);

    if(n_events >= leaf_min_events && n_obs >= leaf_min_obs) {
      cp_min = lc(*it);
      break;
    }

  }

  // reset the counters
  n_events = 0;
  n_obs = 0;

  // beginning from upper values, do the same.

  for(it = lc_sort.end()-1; it >= lc_sort.begin(); --it){

    n_events += y(*it, 1) * weights(*it);
    n_obs += weights(*it);

    if(n_events >= leaf_min_events && n_obs >= leaf_min_obs) {
      cp_max = lc(*it);
      break;
    }


  }

  // if cp_max is < cp_min, neither cut-point is valid.
  // Therefore, we won't split the node because we cant.

  if(cp_max > cp_min){

    Rcout << "cp_max: " << cp_max << std::endl;
    Rcout << "cp_min: " << cp_min << std::endl;

    arma::vec cps_runif = arma::randu<arma::vec>(cp.size());

    for(arma::uword i = 0; i < cp.size(); i++){
      cp(i) = cps_runif(i) * (cp_max - cp_min) + cp_min;
    }

  }


}

// [[Rcpp::export]]
void find_cutpoints(arma::vec& cp,
                    const arma::vec& lc,
                    const arma::mat& y,
                    const arma::uvec weights,
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

  if(lc_size > 50){

    // get unique values in the head and tail
    lc_uni = arma::unique(
      arma::join_cols(lc.head(25),
                      lc.tail(25))
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


        // todo: this could be more efficient
        if(lc[j] <= lc_uni[i]){

          n_events_left += y(j, 1) * weights(j);
          n_obs_left += weights(j);

        } else {

          n_events_right += y(j, 1) * weights(j);
          n_obs_right += weights(j);

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

    if(lc_size > leaf_min_obs * 10){

      arma::uvec rows_sample = arma::linspace<arma::uvec>(0, lc_size-1,
                                                          leaf_min_obs * 10);

      find_cutpoints_ctns(cp,
                          lc(rows_sample),
                          y.rows(rows_sample),
                          weights(rows_sample),
                          leaf_min_obs,
                          leaf_min_events);

    } else {

      find_cutpoints_ctns(cp,
                          lc,
                          y,
                          weights,
                          leaf_min_obs,
                          leaf_min_events);

    }

  }

  //Rcout << cp.t() << std::endl;

}

// [[Rcpp::export]]
double log_rank_test_wtd(arma::mat& y,
                         arma::vec& g,
                         arma::uvec& weights){

  arma::uword n = y.n_rows;
  double Y = arma::sum(weights);
  double Y1 = arma::sum(g % weights);

  arma::uword lwr = 0;
  arma::uword upr = 0;
  arma::uword count = 1; // starts at 1 to mimic size
  arma::uword i;

  for(i=0; i<n; i++){

    if(y(i, 1) == 0){

      upr++;
      count += weights(i);

    } else {

      break;

    }

  }

  double d = 0;
  double d1 = 0;

  for(i = lwr; i <= upr; i++){
    d += y(i, 1) * weights(i);
    d1 += y(i, 1) * weights(i) * g(i);
  }

  double e1 = Y1 * d / Y;
  double e0 = (Y - Y1) * d / Y;
  double o1 = d1;
  double o0 = d - d1;

  double V = (Y-Y1) * Y1 * d * (Y-d) / (pow(Y,2) * (Y-1));

  Y -= count;

  for(i = lwr; i <= upr; i++){
    Y1 -= g(i) * weights(i);
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
    count = weights(upr);


    while( (y(upr, 1) == 0) & (upr < n-1) ){
      upr++;
      count += weights(upr);
    }

    if(upr==n-1){
      if( y(upr, 1) == 0 ){
        break;
      } else {

        d = 0; d1 = 0;

        for(i = lwr; i <= upr; i++){
          d += y(i, 1) * weights(i);
          d1 += y(i, 1) * weights(i) * g(i);
        }

        e1 += (Y1 * d/Y);
        e0 += ((Y-Y1) * d/Y);
        o1 += d1;
        o0 += d-d1;

        V += (Y-Y1) * Y1 * d * (Y-d) / (pow(Y, 2) * (Y-1));
        Y -= count;

        for(i = lwr; i <= upr; i++) Y1 -= g(i) * weights(i);

        break;
      }
    }

    d = 0; d1 = 0;

    for(i = lwr; i <= upr; i++){
      d += y(i, 1) * weights(i);
      d1 += y(i, 1) * weights(i) * g(i);
    }

    e1 += (Y1*d / Y);
    e0 += ((Y-Y1) * d/Y);
    o1 += d1;
    o0 += d - d1;

    V += (Y-Y1) * Y1 * d * (Y-d) / (pow(Y,2) * (Y-1));
    Y -= count;

    for(i = lwr; i <= upr; i++) Y1 -= g(i) * weights(i);

    lwr=upr+1;

    if(Y==1) break;

  }

  return pow(o1-e1,2) / V;

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
double newtraph_cph_iter (const arma::mat& x,
                          const arma::mat& y,
                          const arma::uvec& weights,
                          const arma::vec& beta,
                          arma::vec& XB,
                          arma::vec& R,
                          arma::vec& u,
                          arma::vec& a,
                          arma::vec& a2,
                          arma::mat& imat,
                          arma::mat& cmat,
                          arma::mat& cmat2,
                          const arma::uword& method){

  // know that if you change these to uword it could break the routine
  double wtave, temp, person_time, x_beta, risk, denom2, deadwt, ndead;
  double denom=0, loglik=0;

  arma::uword i, j, k;
  arma::uword nrisk = 0, person = x.n_rows - 1, nvar = x.n_cols;

  u.fill(0);
  a.fill(0);
  a2.fill(0);
  imat.fill(0);
  cmat.fill(0);
  cmat2.fill(0);

  // this loop has a strange break condition to accomodate
  // the restriction that a uvec (uword) cannot be < 0

  bool break_loop = false;

  XB.subvec(0, person) = x * beta;
  R.subvec(0, person) = arma::exp(XB.subvec(0, person)) % weights;

  arma::rowvec x_person(x.n_cols);
  arma::uword weights_person;

  // arma::mat tempmat = x.t() * arma::diagmat(R) * x;
  // Rcout << tempmat << std::endl;


  for( ; ; ){


    // Rcout << "- person: "    << person;
    // Rcout << "; y_status: "  << y(person,1);
    // Rcout << "; y_time: "    << y(person,0);
    // Rcout << "; u: "         << u.t();
    // Rcout << std::endl;

    person_time = y.at(person, 0); // time of event for current person
    ndead  = 0 ; // number of deaths at this time point
    deadwt = 0 ; // sum of weights for the deaths
    denom2 = 0 ; // sum of weighted risks for the deaths

    // walk through this set of tied times
    while(y(person, 0) == person_time){

      nrisk++;

      //x_beta = arma::dot(beta, x.row(person));
      //risk = exp(x_beta) * weights(person);

      x_beta = XB.at(person);
      risk = R.at(person);

      x_person = x.row(person);
      weights_person = weights.at(person);

      if (y.at(person, 1) == 0) {

        denom += risk;

        /* a contains weighted sums of x, cmat sums of squares */

        for (i=0; i<nvar; i++) {

          temp = risk * x_person[i];

          a[i] += temp;

          for (j = 0; j <= i; j++){
            cmat.at(j, i) += temp * x_person[j];
          }

        }

      } else {

        ndead++;

        deadwt += weights_person;
        denom2 += risk;
        loglik += weights_person * x_beta;

        for (i=0; i<nvar; i++) {

          u[i]  += weights_person * x_person[i];
          a2[i] += risk * x_person[i];

          for (j=0; j<=i; j++){
            cmat2.at(j, i) += risk * x_person[i] * x_person[j];
          }

        }

      }

      if(person == 0){
        break_loop = true;
        break;
      }

      person--;

      //if(person_time == y(person, 0)) Rcout << "tie!" << std::endl;

      //Rcout << person << std::endl;

    }

    //Rcout << "imat: " << std::endl << imat << std::endl;

    // we need to add to the main terms
    if (ndead > 0) {

      if (method == 0 || ndead == 1) { // Breslow

        denom  += denom2;
        loglik -= deadwt * log(denom);

        for (i=0; i<nvar; i++) {

          a[i]  += a2[i];
          temp  = a[i] / denom;  // mean
          u[i]  -=  deadwt * temp;

          for (j=0; j<=i; j++) {
            cmat.at(j, i) += cmat2.at(j, i);
            imat.at(j, i) += deadwt * (cmat.at(j, i) - temp * a[j]) / denom;
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
          loglik -= wtave * log(denom);

          for (i=0; i<nvar; i++) {

            a[i] += a2[i] / ndead;
            temp = a[i]  / denom;
            u[i] -= wtave * temp;

            for (j=0; j<=i; j++) {
              cmat.at(j, i) += cmat2.at(j, i) / ndead;
              imat.at(j, i) += wtave * (cmat.at(j, i) - temp * a[j]) / denom;
            }

          }

        }

      }

      a2.fill(0);
      cmat2.fill(0);

    }

    if(break_loop) break;

  }

  //Rcout << cmat << std::endl;

  return(loglik);

}

// [[Rcpp::export]]
arma::mat newtraph_cph(const arma::mat& x,
                       const arma::mat& y,
                       const arma::uvec& weights,
                       const arma::mat& x_transforms,
                       const arma::uword method,
                       const double& eps,
                       const arma::uword iter_max,
                       const bool rescale){

  arma::uword i, j, iter, nvar = x.n_cols;
  double ll_new, ll_best, halving = 0;

  arma::vec beta(nvar);
  arma::vec newbeta(nvar);
  arma::vec u(nvar);
  arma::vec a(nvar);
  arma::vec a2(nvar);

  arma::vec XB(x.n_rows);
  arma::vec R(x.n_rows);

  arma::mat imat(nvar, nvar);
  arma::mat cmat(nvar, nvar);
  arma::mat cmat2(nvar, nvar);

  // do the initial iteration
  ll_best = newtraph_cph_iter(x, y, weights, beta, XB, R, u, a, a2,
                              imat, cmat, cmat2, method);

  //Rcout << ll_best << std::endl;

  // update beta
  cholesky(imat);
  cholesky_solve(imat, u);
  newbeta = beta + u;

  if(iter_max > 1 && !std::isinf(ll_best)){

    for(iter = 1; iter < iter_max; iter++){

      // do the next iteration
      ll_new = newtraph_cph_iter(x, y, weights, newbeta, XB, R, u, a, a2,
                                 imat, cmat, cmat2, method);

      cholesky(imat);

      //Rcout << ll_new << std::endl;

      // check for convergence
      if(abs(1-ll_best/ll_new) < eps){
        break;
      }


      if(ll_new < ll_best){ // it's not converging!

        halving++; // get more aggressive when it doesn't work

        // reduce the magnitude by which newbeta modifies beta
        for (i = 0; i < nvar; i++){
          newbeta[i] = (newbeta[i] + halving * beta[i]) / (halving + 1.0);
        }


      } else { // it's converging!

        halving = 0;
        ll_best = ll_new;
        //Rcout << "log likelihood: " << ll_best << std::endl;

        cholesky_solve(imat, u);

        for (i = 0; i < nvar; i++) {

          beta[i] = newbeta[i];
          newbeta[i] = newbeta[i] +  u[i];

        }

      }

    }

  }

  // invert imat and return to original scale
  cholesky_invert(imat);

  beta = newbeta;

  if(rescale == true){

    for (i=0; i < nvar; i++) {

      beta.at(i) *= x_transforms.at(i, 1);
      u.at(i) /= x_transforms.at(i, 1);
      imat.at(i, i) *= x_transforms.at(i, 1) * x_transforms.at(i, 1);

      for (j = 0; j < i; j++) {

        imat.at(j, i) *= x_transforms.at(i, 1) * x_transforms.at(j, 1);
        imat.at(i, j) = imat.at(j, i);

      }

    }

  }


  for(i = 0; i < nvar; i++){

    if(std::isinf(beta[i])) beta[i] = 0;

    if(std::isinf(imat.at(i, i))) imat.at(i, i) = 1.0;

  }

  arma::vec se = arma::sqrt(imat.diag());

  arma::vec pv(nvar);

  for(i = 0; i < nvar; i++){
    pv[i] = R::pchisq(pow(beta[i]/se[i], 2), 1, false, false);
  }

  arma::mat out(nvar, 3);

  out.col(0) = beta;
  out.col(1) = se;
  out.col(2) = pv;

  return(out);

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
List ostree_fit(arma::mat& x,
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
  // beta coefficients for linear combinations
  arma::mat betas(mtry, max_nodes);
  // column indices to pair with beta coefficients
  arma::umat col_indices(mtry, max_nodes);
  // cutpoints for linear combination splitting
  arma::vec cutpoints(max_nodes);
  // child nodes to progress through the tree
  arma::uvec children_left(max_nodes);

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

  arma::mat x_transforms;
  // vectors and matrices for newton raphson
  arma::vec u(mtry);
  arma::vec a(mtry);
  arma::vec a2(mtry);
  arma::mat imat(mtry, mtry);
  arma::mat cmat(mtry, mtry);
  arma::mat cmat2(mtry, mtry);
  arma::mat cph_fit;


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

  // terms to create new node (nn)
  arma::uword nn_left;

  // label for the new node's name; i.e., node_1, node_2, ...
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

      // if(verbose == true){
      //   Rcout <<
      //     "empty cols_to_sample: " <<
      //       cols_to_sample_01.t() <<
      //         std::endl;
      // }



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

      // if(verbose == true){
      //   Rcout <<
      //     "full  cols_to_sample: " <<
      //       cols_to_sample_01.t() <<
      //         std::endl;
      // }

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

      // if(verbose == true){
      //   Rcout << "cols sampled: " << cols_node.t() << std::endl;
      // }


      x_node = x_inbag(rows_node, cols_node);
      y_node = y_inbag.rows(rows_node);
      weights_node = weights(rows_node);

      if(verbose == true){

        arma::uword n_obs = arma::sum(weights_node);
        arma::uword n_events = arma::sum(y_node.col(1) % weights_node);

        Rcout << "No. of observations in node: " << n_obs << std::endl;
        Rcout << "No. of events in node:       " << n_events << std::endl;

      }

      // Rcout << "x_node: " << std::endl << x_node << std::endl;
      //
      // Rcout << "y_node: " << std::endl << y_node << std::endl;

      //Rcout << "weights: " << weights.t() << std::endl;

      // Rcout << "u: " << u.t() << std::endl;

      x_transforms = x_scale_cph(x_node, weights_node);

      cph_fit = newtraph_cph(x_node,
                             y_node,
                             weights_node,
                             x_transforms,
                             1,      // efron ties
                             1e-4,   // epsilon for convergence
                             3,      // max iterations
                             true);  // rescale coefficients

      // unscale the x matrix
      for(i = 0; i < x_node.n_cols; i++){
        x_node.col(i) /= x_transforms(i,1) + x_transforms(i,0);
      }

      beta = cph_fit.col(0);

      // newtraph_cph_one_iter(x_node, y_node, weights_node,
      //                       u, a, a2, imat, cmat, cmat2, 0);
      //
      // //Rcout << "imat: " << std::endl << imat << std::endl;
      // cholesky(imat);
      // //Rcout << "imat after cholesky: " << std::endl << imat << std::endl;
      // cholesky_solve(imat, u);
      //
      // arma::vec beta(mtry_temp);
      //
      // for(i = 0; i < mtry_temp; i++) beta[i] = u[i];

      if(verbose == true){
        Rcout << "beta: " << beta.t() << std::endl;
      }

      lc = x_node * beta;

      find_cutpoints(cp, lc, y_node, weights_node,
                     leaf_min_obs, leaf_min_events);

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


            lrstat = log_rank_test_wtd(y_node, g, weights_node);

            if(verbose == true){
              Rcout << lrstat << " with cut-point " << cp[i];
            }


            if(lrstat > lrstat_max){
              lrstat_max = lrstat;
              cp_max = cp[i];
              if(verbose == true) Rcout << ", a new max!" << std::endl;
            } else {
              if(verbose == true) Rcout << std::endl;
            }

            if(verbose == true){
              Rcout <<
                "number of g == 0: " <<
                  arma::sum(weights_node) - arma::dot(g, weights_node) <<
                    "; " << arma::sum(y_node.col(1) % weights_node) -
                    arma::sum(y_node.col(1) % weights_node % g) <<
                      " events" << std::endl;
              Rcout << "number of g == 1: " <<
                arma::dot(g, weights_node) <<
                  "; " << arma::sum(y_node.col(1) % weights_node % g) <<
                    " events" << std::endl << std::endl;
            }

            if(verbose == true) Rcout << std::endl;

          }

        }

        nn_left   = nodes_max + 1;
        nodes_max = nodes_max + 2;

        if(verbose == true){
          Rcout << "new nodes: " <<
            nn_left << " and " << nodes_max << std::endl;
        }

        i = 0;

        for(person = rows_node.begin(); person != rows_node.end(); ++person){

          if(lc[i] <= cp_max){

            node_assignments[*person] = nn_left;

          } else {

            node_assignments[*person] = nodes_max;

          }

          i++;

        }

        for(i = 0; i < mtry_temp; i++){
          betas(i, *node) = beta(i);
          col_indices(i, *node) = cols_node(i);
        }

        children_left[*node] = nn_left;
        cutpoints[*node] = cp_max;



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

      arma::uvec temp = arma::regspace<arma::uvec>(0, 1, node_summary.n_rows);

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
        rows_node    = arma::find(node_assignments == i+1);
        y_node       = y_inbag.rows(rows_node);
        weights_node = weights(rows_node);
        leaf         = leaf_surv_small(y_node, weights_node);

        if(verbose == true){
          Rcout << "created new leaf: node_" << i+1 << std::endl;
        }

        leaf_nodes[make_node_name(i+1)] = leaf;

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
      _["mtry"] = mtry
    )
  );




}




// [[Rcpp::export]]
List orsf_fit(arma::mat& x,
              arma::mat& y,
              const int& ntree,
              const arma::uword& mtry = 4,
              const arma::uword& n_cps = 5,
              const arma::uword& leaf_min_events = 5,
              const arma::uword& leaf_min_obs = 10,
              const bool& verbose = false){

  List forest(ntree);

  for(int tree = 0; tree < ntree; tree++){
    forest[tree] = ostree_fit(x,
                              y,
                              mtry,
                              n_cps,
                              leaf_min_events,
                              leaf_min_obs,
                              verbose);
  }

  return(forest);

}

// [[Rcpp::export]]
arma::uvec ostree_pred_leaf(const arma::mat& x_new,
                            const arma::mat& betas,
                            const arma::umat& col_indices,
                            const arma::vec& cut_points,
                            const arma::vec& children_left){

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

            out(k) = children_left(i)+1;

          }

        }

      }

    }

  }

  return(out);

}

// [[Rcpp::export]]
arma::mat ostree_pred_surv(const arma::mat&  x_new,
                           const Rcpp::List& leaf_nodes,
                           const arma::uvec& leaf_preds,
                           const arma::vec&  times){

  // preallocate memory for output
  arma::mat out(x_new.n_rows, times.size());

  arma::uvec leaf_sort = arma::sort_index(leaf_preds);

  //Rcout << leaf_preds(leaf_sort(0)) << std::endl;

  arma::uword person = 0;
  arma::uword person_ref_index;
  arma::uword person_leaf;
  String person_leaf_name;

  arma::uword i, t;

  double surv_estimate;

  do{

    person_ref_index = leaf_sort(person);
    person_leaf = leaf_preds(person_ref_index);

    Rcout << "person: " << person << std::endl;
    // Rcout << "person_ref_index: " << person_ref_index << std::endl;
    // Rcout << "person_leaf: " << person_leaf << std::endl;

    person_leaf_name = make_node_name(person_leaf);

    // got to do it this way to avoid making copy of leaf data
    NumericMatrix leaf_surv_temp = leaf_nodes[person_leaf_name];
    arma::mat leaf_surv(leaf_surv_temp.begin(), leaf_surv_temp.nrow(),
                        leaf_surv_temp.ncol(), false);

    // Rcout << leaf_surv << std::endl;

    i = 0;

    // times must be in ascending order
    // (remember to right a check for this in R API)
    for(t = 0; t < times.size(); t++){

      if(times(t) < leaf_surv(leaf_surv.n_rows - 1, 1)){

        for(; i < leaf_surv.n_rows; i++){
          if (leaf_surv(i, 0) > times(t)){
            if(i == 0)
              surv_estimate = 1;
            else
              surv_estimate = leaf_surv(i-1, 1);
            break;
          } else if (leaf_surv(i, 0) == times(t)){
            surv_estimate = leaf_surv(i, 1);
            break;
          }
        }

      } else {

        // go here if prediction horizon > max time in current leaf.
        double surv_slope =
          (1 - leaf_surv(leaf_surv.n_rows - 1, 1)) /
            (0 - leaf_surv(leaf_surv.n_rows - 1, 0));

        double time_diff = times(t) - leaf_surv(leaf_surv.n_rows - 1, 0);

        surv_estimate =
          leaf_surv(leaf_surv.n_rows - 1, 1) + surv_slope * time_diff;

        if(surv_estimate < 0) surv_estimate = 0;

      }

      out(person_ref_index, t) = surv_estimate;

    }

    person++;

    if(person < x_new.n_rows){

      while(person_leaf == leaf_preds(leaf_sort(person))){

        for(i = 0; i < out.n_cols; i++){
          out(leaf_sort(person), i) = out(person_ref_index, i);
        }

        person++;

        if (person == x_new.n_rows) break;

      }

    }

  } while (person < x_new.n_rows);

  return(out);

}
