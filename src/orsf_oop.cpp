/*-----------------------------------------------------------------------------
 This file is part of aorsf, which is distributed under the MIT license

 You should have received a copy of the MIT License along with aorsf.
 If not, see <http://www.gnu.org/licenses/>.

 Authors:
 - Byron C. Jaeger (http://byronjaeger.com)
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "globals.h"
#include "Data.h"
#include "Tree.h"
#include "Coxph.h"
#include "NodeSplitStats.h"

 // [[Rcpp::depends(RcppArmadillo)]]

 using namespace Rcpp;
 using namespace arma;
 using namespace aorsf;

 // @description sample weights to mimic a bootstrap sample
 // arma::vec bootstrap_sample_testthat(arma::mat& x,
 //                                     arma::mat& y,
 //                                     arma::vec& w) {
 //
 //  Data data = Data(x, y, w);
 //  Data* data_ptr = &data;
 //  return bootstrap_sample(data_ptr);
 //
 // }

 //  Same as x_node_scale, but this can be called from R

 // [[Rcpp::export]]
 List coxph_scale_exported(NumericMatrix& x_,
                           NumericVector& w_){

  mat x_node = mat(x_.begin(), x_.nrow(), x_.ncol(), false);
  vec w_node = vec(w_.begin(), w_.length(), false);
  mat x_transforms = coxph_scale(x_node, w_node);

  return(
   List::create(
    _["x_scaled"] = x_node,
    _["x_transforms"] = x_transforms
   )
  );

 }

 // [[Rcpp::export]]
 List coxph_fit_exported(NumericMatrix& x_,
                         NumericMatrix& y_,
                         NumericVector& w_,
                         int method,
                         double cph_eps,
                         int cph_iter_max){

  mat x_node = mat(x_.begin(), x_.nrow(), x_.ncol(), false);
  mat y_node = mat(y_.begin(), y_.nrow(), y_.ncol(), false);
  vec w_node = vec(w_.begin(), w_.length(), false);

  uword cph_iter_max_ = cph_iter_max;


  vec beta = coxph_fit(x_node,
                       y_node,
                       w_node,
                       method,
                       cph_eps,
                       cph_iter_max_,
                       'A');

  return(
   List::create(
    _["beta"] = beta
   )
  );

 }

 // [[Rcpp::export]]
 arma::vec node_find_cps(arma::mat& y_node,
                         arma::vec& w_node,
                         arma::vec& XB,
                         arma::uword n_split,
                         double leaf_min_events,
                         double leaf_min_obs){


  vec
   // unsafe_cols point to cols in y_node.
   y_time = y_node.unsafe_col(0),
   y_status = y_node.unsafe_col(1),
   cp_placeholder(n_split),
   cp_lwr(1),
   cp_mid,
   cp_upr(1);

  // sort XB to iterate over the sorted indices
  uvec XB_sorted = sort_index(XB, "ascend");

  cp_placeholder.fill(R_PosInf);

  uword i, j, k;

  uvec::iterator
   XB_iter,
   XB_iter_lwr,
   XB_iter_upr;

  double n_events = 0, n_risk = 0;

  if(VERBOSITY > 0){
   Rcout << "----- finding a lower bound for cut-points -----" << std::endl;
  }

  // stop at end-1 b/c we access XB_iter+1 in XB_sorted
  for(XB_iter = XB_sorted.begin(); XB_iter < XB_sorted.end()-1; ++XB_iter){

   n_events += y_status[*XB_iter] * w_node[*XB_iter];
   n_risk += w_node[*XB_iter];


   if(VERBOSITY > 1){
    Rcout << "current XB"<< XB(*XB_iter)  << " ---- ";
    Rcout << "next XB"<< XB(*(XB_iter+1)) << " ---- ";
    Rcout << "N events" << n_events       << " ---- ";
    Rcout << "N risk" << n_risk           << std::endl;
   }

   // If we want to make the current value of XB a cut-point, we need
   // to make sure the next value of XB isn't equal to this current value.
   // Otherwise, we will have the same value of XB in both groups!

   if(XB[*XB_iter] != XB[*(XB_iter+1)]){

    if(VERBOSITY > 1){
     Rcout << "********* New cut-point here ********" << std::endl;
    }

    if( n_events >= leaf_min_events &&
        n_risk   >= leaf_min_obs) {

     if(VERBOSITY > 1){
      Rcout << std::endl;
      Rcout << "lower cutpoint: "         << XB(*XB_iter) << std::endl;
      Rcout << " - n_events, left node: " << n_events << std::endl;
      Rcout << " - n_risk, left node:   " << n_risk   << std::endl;
      Rcout << std::endl;
     }

     break;

    }

   }

  }

  if(XB_iter == XB_sorted.end()-1) {

   if(VERBOSITY > 1){
    Rcout << "Could not find a valid lower cut-point" << std::endl;
   }

   return(cp_placeholder);

  }

  XB_iter_lwr = XB_iter;
  cp_lwr[0] = XB[*XB_iter];

  // set j to be the number of steps we have taken forward in XB
  j = XB_iter - XB_sorted.begin();

  // reset before finding the upper limit
  n_events=0, n_risk=0;

  // stop at beginning+1 b/c we access XB_iter-1 in XB_sorted
  for(XB_iter = XB_sorted.end()-1; XB_iter >= XB_sorted.begin()+1; --XB_iter){

   n_events += y_status[*XB_iter] * w_node[*XB_iter];
   n_risk   += w_node[*XB_iter];

   if(VERBOSITY > 1){
    Rcout << XB(*XB_iter)     << " ---- ";
    Rcout << XB(*(XB_iter-1)) << " ---- ";
    Rcout << n_events     << " ---- ";
    Rcout << n_risk       << std::endl;
   }

   if(XB(*XB_iter) != XB(*(XB_iter-1))){

    if(VERBOSITY > 1){
     Rcout << "********* New cut-point here ********" << std::endl;
    }

    if( n_events >= leaf_min_events &&
        n_risk   >= leaf_min_obs ) {

     // the upper cutpoint needs to be one step below the current
     // XB_iter value, because we use x <= cp to determine whether a
     // value x goes to the left node versus the right node. So,
     // if XB_iter currently points to 3, and the next value down is 2,
     // then we want to say the cut-point is 2 because then all
     // values <= 2 will go left, and 3 will go right. This matters
     // when 3 is the highest value in the vector.

     --XB_iter;

     if(VERBOSITY > 1){
      Rcout << std::endl;
      Rcout << "upper cutpoint: " << XB(*XB_iter) << std::endl;
      Rcout << " - n_events, right node: " << n_events    << std::endl;
      Rcout << " - n_risk, right node:   " << n_risk      << std::endl;
     }

     break;

    }

   }

  }

  // k = n steps from beginning of sorted XB to current XB_iter
  k = XB_iter + 1 - XB_sorted.begin();

  if(VERBOSITY > 1){
   Rcout << "N steps from beginning to first cp: " << j << std::endl;
   Rcout << "N steps from beginning to last cp: " << k << std::endl;
   Rcout << "n potential cutpoints: " << k-j << std::endl;
  }

  if(j > k){

   if(VERBOSITY > 1) {
    Rcout << "Could not find valid cut-points" << std::endl;
   }

   return(cp_placeholder);

  }

  XB_iter_upr = XB_iter;
  cp_upr[0] = XB[*XB_iter];

  uvec tmp = XB_sorted(*XB_iter_lwr, *XB_iter_upr);
  Rcout << XB(tmp) << std::endl;

  if(cp_lwr[0] != cp_upr[0]){

  }

  return(join_vert(cp_lwr, cp_mid, cp_upr));


 }

 // [[Rcpp::export]]
 List lrt_multi_exported(NumericMatrix& y_,
                         NumericVector& w_,
                         NumericVector& XB_,
                         int n_split_,
                         double split_min_stat,
                         double leaf_min_events,
                         double leaf_min_obs){

  mat y_node = mat(y_.begin(), y_.nrow(), y_.ncol(), false);
  vec w_node = vec(w_.begin(), w_.length(), false);
  vec XB = vec(XB_.begin(), XB_.length(), false);

  uword n_split = n_split_;

  List out = lrt_multi(y_node,
                       w_node,
                       XB,
                       n_split,
                       split_min_stat,
                       leaf_min_events,
                       leaf_min_obs);

  return(out);

 }


 // [[Rcpp::export]]
 List orsf_cpp(arma::mat& x,
               arma::mat& y,
               arma::vec& w,
               int vi = 0,
               int sr = 1,
               int pt = 1,
               bool oobag_pred = true){


  int mtry = 2;
  int leaf_min_obs = DEFAULT_LEAF_MIN_OBS_SURVIVAL;
  // int split_min_obs = DEFAULT_SPLIT_MIN_OBS;
  // int split_min_stat = DEFAULT_SPLIT_MIN_STAT;
  // int max_retry = DEFAULT_MAX_RETRY;
  // int n_split = DEFAULT_N_SPLIT;
  // int oobag_eval_every = 0;
  // int seed = 0;

  VariableImportance variable_importance = static_cast<VariableImportance>(vi);
  // SplitRule split_rule = static_cast<SplitRule>(sr);
  // PredType pred_type = static_cast<PredType>(pt);

  // if( variable_importance == VI_NONE )
  //  Rcout << variable_importance << std::endl;

  Data data = Data(x, y, w);

  if(VERBOSITY > 0){
   Rcout << "------------ dimensions ------------"   << std::endl;
   Rcout << "N obs total: "     << data.get_n_rows() << std::endl;
   Rcout << "N columns total: " << data.get_n_cols() << std::endl;
   Rcout << "------------------------------------";
   Rcout << std::endl << std::endl;
  }

  Rcpp::List result;

  Data* data_ptr = &data;

  Tree tree(data_ptr, leaf_min_obs, mtry);

  tree.grow(oobag_pred);

  result.push_back(tree.rows_oobag, "rows_oobag");
  result.push_back(tree.coef, "coef");
  result.push_back(tree.coef_indices, "coef_indices");
  result.push_back(tree.cutpoint, "cutpoint");
  result.push_back(tree.next_left_node, "next_left_node");
  result.push_back(tree.leaf_values, "leaf_values");
  result.push_back(tree.pred_oobag, "pred_oobag");
  result.push_back(tree.leaf_indices, "leaf_indices");


  return(result);


 }
