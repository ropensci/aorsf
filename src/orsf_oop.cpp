/*-----------------------------------------------------------------------------
 This file is part of aorsf, which is distributed under the MIT license

 You should have received a copy of the MIT License along with aorsf.
 If not, see <http://www.gnu.org/licenses/>.

 Authors:
 - Byron C. Jaeger (http://byronjaeger.com)
 - test
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
 List coxph_scale_exported(arma::vec& x_node,
                           arma::vec& w_node){

  mat x_transforms = coxph_scale(x_node, w_node);

  return(
   List::create(
    _["x_scaled"] = x_node,
    _["x_transforms"] = x_transforms
   )
  );

 }

 // [[Rcpp::export]]
 List coxph_fit_exported(arma::mat& x_node,
                         arma::mat& y_node,
                         arma::vec& w_node,
                         int method,
                         double cph_eps,
                         arma::uword cph_iter_max){

  arma::vec vi_numer(x_node.n_cols);
  arma::uvec vi_denom(x_node.n_cols);
  arma::uvec cols_node=regspace<uvec>(0, x_node.n_cols-1);

  vec beta = coxph_fit(x_node,
                       y_node,
                       w_node,
                       cols_node,
                       true,
                       method,
                       cph_eps,
                       cph_iter_max,
                       0.10,
                       VI_ANOVA,
                       vi_numer,
                       vi_denom);

  return(
   List::create(
    _["beta"] = beta,
    _["vi_numer"] = vi_numer,
    _["vi_denom"] = vi_denom
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
                               arma::uvec& XB_sorted,
                               arma::uword start,
                               arma::uword stop,
                               double value){

  node_fill_group(group, XB_sorted, start, stop, value);

 }

 // [[Rcpp::plugins("cpp17")]]
 // [[Rcpp::export]]
 List orsf_cpp(arma::mat& x,
               arma::mat& y,
               arma::vec& w,
               Rcpp::IntegerVector& tree_seeds,
               Rcpp::Function& lincomb_R_function,
               Rcpp::Function& oobag_R_function,
               arma::uword n_tree,
               arma::uword mtry,
               arma::uword vi_type_R,
               double leaf_min_events,
               double leaf_min_obs,
               arma::uword split_rule_R,
               double split_min_events,
               double split_min_obs,
               double split_min_stat,
               arma::uword split_max_cuts,
               arma::uword split_max_retry,
               arma::uword lincomb_type_R,
               double lincomb_eps,
               arma::uword lincomb_iter_max,
               bool lincomb_scale,
               double lincomb_alpha,
               arma::uword lincomb_df_target,
               arma::uword lincomb_ties_method,
               bool pred_mode,
               arma::uword pred_type_R,
               double pred_horizon,
               bool oobag_pred,
               arma::uword oobag_eval_every){

  // int mtry = 2;
  // int leaf_min_obs = DEFAULT_LEAF_MIN_OBS_SURVIVAL;
  // VariableImportance variable_importance = static_cast<VariableImportance>(vi);
  // SplitRule split_rule = static_cast<SplitRule>(sr);
  // PredType pred_type = static_cast<PredType>(pt);

  // VariableImportance vi_type = static_cast<VariableImportance>(vi_type_R);
  // SplitRule split_rule = static_cast<SplitRule>(split_rule_r) ;
  // LinearCombo lincomb_type = static_cast<LinearCombo>(lincomb_type_R);
  // PredType pred_type = static_cast<PredType>(pred_type_R);

  Rcpp::List result, forest_out;

  std::unique_ptr<Forest> forest { };
  std::unique_ptr<Data> data { };

  data = std::make_unique<Data>(x, y, w);

  VariableImportance vi_type = (VariableImportance) vi_type_R;
  SplitRule split_rule = (SplitRule) split_rule_R;
  LinearCombo lincomb_type = (LinearCombo) lincomb_type_R;
  PredType pred_type = (PredType) pred_type_R;

  forest = std::make_unique<Forest>();

  forest->init(std::move(data),
               tree_seeds,
               n_tree,
               mtry,
               vi_type,
               leaf_min_events,
               leaf_min_obs,
               split_rule,
               split_min_events,
               split_min_obs,
               split_min_stat,
               split_max_cuts,
               split_max_retry,
               lincomb_type,
               lincomb_eps,
               lincomb_iter_max,
               lincomb_scale,
               lincomb_alpha,
               lincomb_df_target,
               lincomb_ties_method,
               pred_type,
               pred_mode,
               pred_horizon,
               oobag_pred,
               oobag_eval_every);

  forest->plant();

  forest->grow(lincomb_R_function);

  forest_out.push_back(forest->get_coef_indices(), "coef_indices");
  forest_out.push_back(forest->get_leaf_pred_horizon(), "leaf_pred_horizon");
  forest_out.push_back(forest->get_leaf_pred_surv(), "leaf_pred_surv");
  forest_out.push_back(forest->get_leaf_pred_chf(), "leaf_pred_chf");


  result.push_back(forest_out, "forest");

  return(result);

 }
