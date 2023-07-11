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
