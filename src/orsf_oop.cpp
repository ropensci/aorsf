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
#include "Forest.h"
#include "Coxph.h"
#include "NodeSplitStats.h"

#include <utility>
#include <memory>
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
 List node_find_cps_exported(arma::mat& y_node,
                             arma::vec& w_node,
                             arma::vec& XB,
                             double leaf_min_events,
                             double leaf_min_obs){

  // sort XB to iterate over the sorted indices
  uvec XB_sorted = sort_index(XB, "ascend");

  uvec cp_index = node_find_cps(y_node,
                                w_node,
                                XB,
                                XB_sorted,
                                leaf_min_events,
                                leaf_min_obs);

  // vec group(y_node.n_rows, fill::zeros);
  // uvec::iterator XB_iter;
  // uword j = 0;
  // XB_iter = XB_sorted.begin();
  // while(j <= cp_index(0)){
  //  group(*XB_iter) = 1;
  //  ++XB_iter;
  //  ++j;
  // }



  return(
   List::create(
    _["cp_index"] = cp_index,
    _["XB_sorted"] = XB_sorted
   )
  );


 }

 // [[Rcpp::export]]
 double node_compute_lrt_exported(arma::mat& y_node,
                                  arma::vec& w_node,
                                  arma::vec& group){

  double out = node_compute_lrt(y_node, w_node, group);

  return(out);

 }

 // [[Rcpp::export]]
 void node_fill_group_exported(arma::vec& group,
                               const arma::uvec& XB_sorted,
                               const arma::uword start,
                               const arma::uword stop,
                               const double value){

  node_fill_group(group, XB_sorted, start, stop, value);

 }

 // valid columns are non-constant in the rows where events occurred
 // [[Rcpp::export]]
 arma::uvec which_cols_valid_exported(const arma::mat& y_inbag,
                                      const arma::mat& x_inbag,
                                      arma::uvec& rows_node,
                                      const arma::uword mtry){

  uvec result(x_inbag.n_cols, arma::fill::zeros);

  // j moves along columns, i along rows, and iit along
  uword j;
  uvec::iterator iit;

  double temp1;//, temp2;

  for(j = 0; j < result.size(); j++){

   temp1 = R_PosInf;

   for(iit = rows_node.begin(); iit != rows_node.end(); ++iit){

    if(y_inbag.at(*iit, 1) == 1){

     if (temp1 < R_PosInf){

      if(x_inbag.at(*iit, j) != temp1){

       result[j] = 1;
       break;

      }

     } else {

      temp1 = x_inbag.at(*iit, j);

     }

    }

   }

  }

  return(result);

 }


 // deprecated, need to drop this
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


 // [[Rcpp::plugins("cpp17")]]
 // [[Rcpp::export]]

 List orsf_cpp(arma::mat& x,
               arma::mat& y,
               arma::vec& w,
               int n_tree,
               Rcpp::Function f_beta,
               Rcpp::Function f_oobag_eval,
               Rcpp::IntegerVector& tree_seeds,
               Rcpp::List& tree_params){


  // int mtry = 2;
  // int leaf_min_obs = DEFAULT_LEAF_MIN_OBS_SURVIVAL;
  // VariableImportance variable_importance = static_cast<VariableImportance>(vi);
  // SplitRule split_rule = static_cast<SplitRule>(sr);
  // PredType pred_type = static_cast<PredType>(pt);


  std::unique_ptr<Forest> forest { };
  std::unique_ptr<Data> data { };

  data = std::make_unique<Data>(x, y, w);

  forest = std::make_unique<Forest>();

  forest->init(std::move(data),
               n_tree,
               tree_seeds,
               tree_params);

  Rcpp::List result;

  result.push_back(
   forest->get_bootstrap_select_probs(),
   "bootstrap_select_probs"
  );


  return(result);


 }
