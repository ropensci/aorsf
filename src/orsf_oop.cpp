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
#include "ForestSurvival.h"
#include "Coxph.h"
#include "NodeSplitStats.h"

#include <utility>
#include <memory>

// [[Rcpp::depends(RcppArmadillo)]]

 using namespace Rcpp;
 using namespace arma;
 using namespace aorsf;

 // [[Rcpp::export]]
 List coxph_fit_exported(arma::mat& x_node,
                         arma::mat& y_node,
                         arma::vec& w_node,
                         int method,
                         double cph_eps,
                         arma::uword cph_iter_max){

  arma::uvec cols_node=regspace<uvec>(0, x_node.n_cols-1);

  arma::mat out = coxph_fit(x_node,
                            y_node,
                            w_node,
                            true,
                            method,
                            cph_eps,
                            cph_iter_max);

  List result;
  result.push_back(out.col(0), "beta");
  result.push_back(out.col(1), "var");

  return(result);

 }

 // [[Rcpp::export]]
 double compute_cstat_exported_vec(
   arma::mat& y,
   arma::vec& w,
   arma::vec& p,
   bool pred_is_risklike
 ){ return compute_cstat(y, w, p, pred_is_risklike); }

 // [[Rcpp::export]]
 double compute_cstat_exported_uvec(
   arma::mat& y,
   arma::vec& w,
   arma::uvec& g,
   bool pred_is_risklike
 ){ return compute_cstat(y, w, g, pred_is_risklike); }

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
               arma::uword tree_type_R,
               Rcpp::IntegerVector& tree_seeds,
               Rcpp::List loaded_forest,
               Rcpp::RObject lincomb_R_function,
               Rcpp::RObject oobag_R_function,
               arma::uword n_tree,
               arma::uword mtry,
               arma::uword vi_type_R,
               double vi_max_pvalue,
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
               arma::vec pred_horizon,
               bool oobag,
               arma::uword oobag_eval_every,
               unsigned int n_thread){

  List result, forest_out;

  std::unique_ptr<Forest> forest { };
  std::unique_ptr<Data> data { };

  data = std::make_unique<Data>(x, y, w);

  // re-cast integer inputs from R into enumerations
  // see globals.h for definitions.
  TreeType tree_type = (TreeType) tree_type_R;
  VariableImportance vi_type = (VariableImportance) vi_type_R;
  SplitRule split_rule = (SplitRule) split_rule_R;
  LinearCombo lincomb_type = (LinearCombo) lincomb_type_R;
  PredType pred_type = (PredType) pred_type_R;

  // R functions cannot be called from multiple threads
  if(lincomb_type == LC_R_FUNCTION || lincomb_type == LC_GLMNET){
   n_thread = 1;
  }

  // usually need to set n_thread to 1 if oobag pred is monitored
  if(oobag_eval_every < n_tree){
   // specifically if this isn't true we need to go single thread
   if(n_tree/oobag_eval_every != n_thread){
    n_thread = 1;
   }
  }

  if(tree_type == TREE_SURVIVAL){

   vec unique_event_times = find_unique_event_times(y);

   forest = std::make_unique<ForestSurvival>(leaf_min_events,
                                             split_min_events,
                                             pred_horizon,
                                             unique_event_times);

  } else {

   Rcpp::stop("only survival trees are currently implemented");

  }

  forest->init(std::move(data),
               tree_seeds,
               n_tree,
               mtry,
               vi_type,
               vi_max_pvalue,
               leaf_min_obs,
               split_rule,
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
               lincomb_R_function,
               pred_type,
               pred_mode,
               oobag,
               oobag_eval_every,
               oobag_R_function,
               n_thread);

   // Load forest object if in prediction mode
  if(pred_mode){

   std::vector<std::vector<double>> cutpoint     = loaded_forest["cutpoint"];
   std::vector<std::vector<uword>>  child_left   = loaded_forest["child_left"];
   std::vector<std::vector<vec>>    coef_values  = loaded_forest["coef_values"];
   std::vector<std::vector<uvec>>   coef_indices = loaded_forest["coef_indices"];
   std::vector<std::vector<double>> leaf_summary = loaded_forest["leaf_summary"];

   if(tree_type == TREE_SURVIVAL){

    std::vector<std::vector<vec>> leaf_pred_indx = loaded_forest["leaf_pred_indx"];
    std::vector<std::vector<vec>> leaf_pred_prob = loaded_forest["leaf_pred_prob"];
    std::vector<std::vector<vec>> leaf_pred_chaz = loaded_forest["leaf_pred_chaz"];

    auto& temp = dynamic_cast<ForestSurvival&>(*forest);
    temp.load(n_tree, cutpoint, child_left, coef_values, coef_indices,
              leaf_pred_indx, leaf_pred_prob, leaf_pred_chaz, leaf_summary);

   }

   arma::mat pred_mat = forest->predict(oobag);

   result.push_back(pred_mat, "predictions");

  } else {

   // initialize the trees
   forest->plant();

   // grow the trees
   forest->grow();

   // compute out-of-bag predictions if needed
   if(oobag){

    mat pred_oobag = forest->predict(oobag);

    result.push_back(pred_oobag, "pred_oobag");

    if(oobag_eval_every == n_tree){
     forest->compute_prediction_accuracy(y, w, 0, pred_oobag);
    }

    List eval_oobag_out = List::create(
     _["stat_values"] = forest->get_oobag_eval(),
     _["stat_type"] = 1
    );

    result.push_back(eval_oobag_out, "eval_oobag");
   }


   forest_out.push_back(forest->get_rows_oobag(), "rows_oobag");
   forest_out.push_back(forest->get_cutpoint(), "cutpoint");
   forest_out.push_back(forest->get_child_left(), "child_left");
   forest_out.push_back(forest->get_coef_indices(), "coef_indices");
   forest_out.push_back(forest->get_coef_values(), "coef_values");
   forest_out.push_back(forest->get_leaf_summary(), "leaf_summary");

   if(tree_type == TREE_SURVIVAL){
    auto& temp = dynamic_cast<ForestSurvival&>(*forest);
    forest_out.push_back(temp.get_leaf_pred_indx(), "leaf_pred_indx");
    forest_out.push_back(temp.get_leaf_pred_prob(), "leaf_pred_prob");
    forest_out.push_back(temp.get_leaf_pred_chaz(), "leaf_pred_chaz");
    result.push_back(forest->get_unique_event_times(), "unique_event_times");
    result.push_back(pred_horizon, "pred_horizon");
   }

   result.push_back(forest_out, "forest");

   vec vi_output;

   if(vi_type != VI_NONE){

    if(vi_type == VI_ANOVA){
     vi_output = forest->get_vi_numer() / forest->get_vi_denom();
    } else {
     vi_output = forest->get_vi_numer() / n_tree;
    }

   }

   result.push_back(vi_output, "importance");

  }

  return(result);

 }
