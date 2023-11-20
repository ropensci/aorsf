/*-----------------------------------------------------------------------------
 This file is part of aorsf, which is distributed under the MIT license

 You should have received a copy of the MIT License along with aorsf.
 If not, see <http://www.gnu.org/licenses/>.

 Authors:
 - Byron C. Jaeger (http://byronjaeger.com)
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include <vector>
#include <memory>
#include <utility>

#include "globals.h"
#include "Data.h"
#include "Tree.h"
#include "Forest.h"
#include "ForestSurvival.h"
#include "ForestClassification.h"
#include "ForestRegression.h"
#include "Coxph.h"
#include "utility.h"

// [[Rcpp::depends(RcppArmadillo)]]

 using namespace Rcpp;
 using namespace arma;
 using namespace aorsf;

 // [[Rcpp::export]]
 List coxph_fit_exported(arma::mat& x_node,
                         arma::mat& y_node,
                         arma::vec& w_node,
                         int method,
                         double epsilon,
                         arma::uword iter_max){

  arma::uvec cols_node=regspace<uvec>(0, x_node.n_cols-1);

  arma::mat out = coxph_fit(x_node,
                            y_node,
                            w_node,
                            true,
                            method,
                            epsilon,
                            iter_max);

  List result;
  result.push_back(out.col(0), "beta");
  result.push_back(out.col(1), "pvalues");

  return(result);

 }

 // [[Rcpp::export]]
 arma::mat linreg_fit_exported(arma::mat& x_node,
                               arma::mat& y_node,
                               arma::vec& w_node,
                               bool do_scale,
                               double epsilon,
                               arma::uword iter_max){
  return(
   linreg_fit(x_node, y_node, w_node, do_scale, epsilon, iter_max)
  );
 }

 // [[Rcpp::export]]
 arma::mat logreg_fit_exported(arma::mat& x_node,
                               arma::mat& y_node,
                               arma::vec& w_node,
                               bool do_scale,
                               double epsilon,
                               arma::uword iter_max){
  return(
   logreg_fit(x_node, y_node, w_node, do_scale, epsilon, iter_max)
  );
 }

 // [[Rcpp::export]]
 double compute_cstat_surv_exported_vec(
   arma::mat& y,
   arma::vec& w,
   arma::vec& p,
   bool pred_is_risklike
 ){ return compute_cstat_surv(y, w, p, pred_is_risklike); }

 // [[Rcpp::export]]
 double compute_cstat_surv_exported_uvec(
   arma::mat& y,
   arma::vec& w,
   arma::uvec& g,
   bool pred_is_risklike
 ){ return compute_cstat_surv(y, w, g, pred_is_risklike); }

 // [[Rcpp::export]]
 double compute_cstat_clsf_exported(
   arma::vec& y,
   arma::vec& w,
   arma::vec& p
 ){ return compute_cstat_clsf(y, w, p); }

 // [[Rcpp::export]]
 double compute_logrank_exported(
   arma::mat& y,
   arma::vec& w,
   arma::uvec& g
 ){ return compute_logrank(y, w, g); }

 // [[Rcpp::export]]
 double compute_gini_exported(
   arma::mat& y,
   arma::vec& w,
   arma::uvec& g
 ){ return compute_gini(y, w, g); }

 // [[Rcpp::export]]
 arma::vec compute_pred_prob_exported(
   arma::mat& y,
   arma::vec& w
 ){ return compute_pred_prob(y, w); }

 // [[Rcpp::export]]
 double compute_var_reduction_exported(arma::vec& y_node,
                                       arma::vec& w_node,
                                       arma::uvec& g_node){

  return(compute_var_reduction(y_node, w_node, g_node));

 }


 // [[Rcpp::export]]
 bool is_col_splittable_exported(arma::mat& x,
                                 arma::mat& y,
                                 arma::uvec& r,
                                 arma::uword j){

  TreeSurvival tree;

  tree.set_x_inbag(x);
  tree.set_y_inbag(y);
  tree.set_rows_node(r);

  return(tree.is_col_splittable(j));

 }

 // [[Rcpp::export]]
 List find_cuts_survival_exported(arma::mat& y,
                                  arma::vec& w,
                                  arma::vec& lincomb,
                                  double leaf_min_events,
                                  double leaf_min_obs,
                                  int split_rule_R){


  TreeSurvival tree;
  SplitRule split_rule = (SplitRule) split_rule_R;
  arma::uvec lincomb_sort = sort_index(lincomb);

  tree.set_y_node(y);
  tree.set_w_node(w);
  tree.set_lincomb(lincomb);
  tree.set_lincomb_sort(lincomb_sort);
  tree.set_leaf_min_obs(leaf_min_obs);
  tree.set_leaf_min_events(leaf_min_events);
  tree.set_seed(329);
  tree.set_split_max_cuts(5);
  tree.set_split_rule(split_rule);

  tree.find_all_cuts();
  tree.sample_cuts();

  double best_cut = tree.find_best_cut();

  List result;

  result.push_back(tree.get_cuts_all(), "cuts_all");
  result.push_back(tree.get_cuts_sampled(), "cuts_sampled");
  result.push_back(best_cut, "best_cut");

  return(result);

 }

 // [[Rcpp::export]]
 List sprout_node_survival_exported(arma::mat& y,
                                    arma::vec& w){

  TreeSurvival tree;

  arma::uword node_id = 0;
  arma::uword leaf_size = 1;

  tree.unique_event_times = new vec(find_unique_event_times(y));
  tree.set_y_node(y);
  tree.set_w_node(w);

  tree.resize_leaves(leaf_size);
  tree.sprout_leaf(node_id);

  List result;

  result.push_back(tree.get_leaf_pred_indx(), "indx");
  result.push_back(tree.get_leaf_pred_prob(), "prob");
  result.push_back(tree.get_leaf_pred_chaz(), "chaz");
  result.push_back(tree.get_leaf_summary(),   "mort");

  delete tree.unique_event_times;

  return(result);

 }

 // [[Rcpp::export]]
 arma::uvec find_rows_inbag_exported(arma::uvec rows_oobag,
                                     arma::uword n_obs){

  TreeSurvival tree;

  tree.set_rows_oobag(rows_oobag);
  tree.find_rows_inbag(n_obs);

  return(tree.get_rows_inbag());

 }

 // [[Rcpp::export]]
 arma::vec x_submat_mult_beta_exported(arma::mat& x,
                                       arma::mat& y,
                                       arma::vec& w,
                                       arma::uvec& x_rows,
                                       arma::uvec& x_cols,
                                       arma::vec& beta){

  std::unique_ptr<Data> data = std::make_unique<Data>(x, y, w);

  vec out = data->x_submat_mult_beta(x_rows, x_cols, beta);

  return(out);

 }

 // [[Rcpp::export]]
 arma::vec x_submat_mult_beta_pd_exported(arma::mat& x,
                                          arma::mat& y,
                                          arma::vec& w,
                                          arma::uvec& x_rows,
                                          arma::uvec& x_cols,
                                          arma::vec& beta,
                                          arma::vec& pd_x_vals,
                                          arma::uvec& pd_x_cols){

  std::unique_ptr<Data> data = std::make_unique<Data>(x, y, w);

  vec out = data->x_submat_mult_beta(x_rows, x_cols, beta,
                                     pd_x_vals, pd_x_cols);

  return(out);

 }

 // [[Rcpp::export]]
 List scale_x_exported(arma::mat& x,
                       arma::vec& w){

  mat transforms = scale_x(x, w);

  List result;
  result.push_back(x, "scaled_x");
  result.push_back(transforms, "transforms");

  return(result);

 }

 // [[Rcpp::export]]
 List cph_scale(arma::mat& x,
                arma::vec& w){

  // set aside memory for outputs
  // first column holds the mean values
  // second column holds the scale values

  uword n_vars = x.n_cols;

  mat x_transforms(n_vars, 2, fill::zeros);
  vec means  = x_transforms.unsafe_col(0);   // Reference to column 1
  vec scales = x_transforms.unsafe_col(1);   // Reference to column 2

  double w_sum = sum(w);

  for(uword i = 0; i < n_vars; i++) {

   means.at(i) = sum( w % x.col(i) ) / w_sum;

   x.col(i) -= means.at(i);

   scales.at(i) = sum(w % abs(x.col(i)));

   if(scales(i) > 0)
    scales.at(i) = w_sum / scales.at(i);
   else
    scales.at(i) = 1.0; // rare case of constant covariate;

   x.col(i) *= scales.at(i);

  }

  List result;
  result.push_back(x, "x_scaled");
  result.push_back(x_transforms, "x_transforms");
  return result;

 }

 // [[Rcpp::export]]
 arma::mat expand_y_clsf(arma::vec& y,
                         arma::uword n_class){

  arma::mat out(y.n_rows, n_class, arma::fill::zeros);

  for(arma::uword i = 0; i < y.n_rows; ++i){
   out.at(i, y[i]) = 1;
  }

  return(out);

 }

// [[Rcpp::export]]
double compute_mse_exported(arma::vec& y,
                            arma::vec& w,
                            arma::vec& p){

 return(compute_mse(y, w, p));

}


 // [[Rcpp::export]]
 List orsf_cpp(arma::mat&               x,
               arma::mat&               y,
               arma::vec&               w,
               arma::uword              tree_type_R,
               Rcpp::IntegerVector&     tree_seeds,
               Rcpp::List&              loaded_forest,
               Rcpp::RObject            lincomb_R_function,
               Rcpp::RObject            oobag_R_function,
               arma::uword              n_tree,
               arma::uword              mtry,
               bool                     sample_with_replacement,
               double                   sample_fraction,
               arma::uword              vi_type_R,
               double                   vi_max_pvalue,
               double                   leaf_min_events,
               double                   leaf_min_obs,
               arma::uword              split_rule_R,
               double                   split_min_events,
               double                   split_min_obs,
               double                   split_min_stat,
               arma::uword              split_max_cuts,
               arma::uword              split_max_retry,
               arma::uword              lincomb_type_R,
               double                   lincomb_eps,
               arma::uword              lincomb_iter_max,
               bool                     lincomb_scale,
               double                   lincomb_alpha,
               arma::uword              lincomb_df_target,
               arma::uword              lincomb_ties_method,
               bool                     pred_mode,
               arma::uword              pred_type_R,
               arma::vec                pred_horizon,
               bool                     pred_aggregate,
               bool                     oobag,
               arma::uword              oobag_eval_type_R,
               arma::uword              oobag_eval_every,
               int                      pd_type_R,
               std::vector<arma::mat>&  pd_x_vals,
               std::vector<arma::uvec>& pd_x_cols,
               arma::vec&               pd_probs,
               unsigned int             n_thread,
               bool                     write_forest,
               bool                     run_forest,
               int                      verbosity){

  // re-cast integer inputs from R into enumerations
  // see globals.h for definitions.
  VariableImportance vi_type = (VariableImportance) vi_type_R;
  SplitRule split_rule = (SplitRule) split_rule_R;
  LinearCombo lincomb_type = (LinearCombo) lincomb_type_R;
  PredType pred_type = (PredType) pred_type_R;
  EvalType oobag_eval_type = (EvalType) oobag_eval_type_R;
  PartialDepType pd_type = (PartialDepType) pd_type_R;
  TreeType tree_type = (TreeType) tree_type_R;

  List result;

  std::unique_ptr<Forest> forest { };
  std::unique_ptr<Data> data { };

  data = std::make_unique<Data>(x, y, w);

  uword n_obs = data->get_n_rows();

  if(n_thread == 0){
   n_thread = std::thread::hardware_concurrency();
  }

  // R functions cannot be called from multiple threads
  if(lincomb_type    == LC_R_FUNCTION  ||
     lincomb_type    == LC_GLMNET      ||
     oobag_eval_type == EVAL_R_FUNCTION){
   n_thread = 1;
  }

  // usually need to set n_thread to 1 if oobag pred is monitored
  if(oobag_eval_every < n_tree){
   // specifically if this isn't true we need to go single thread
   if(n_tree/oobag_eval_every != n_thread){
    n_thread = 1;
   }
  }

  switch(tree_type){

  case TREE_SURVIVAL:

   forest = std::make_unique<ForestSurvival>(leaf_min_events,
                                             split_min_events,
                                             pred_horizon);

   if(verbosity > 3){
    Rcout << "initializing survival forest" << std::endl;
    Rcout << "  -- leaf_min_events: " << leaf_min_events << std::endl;
    Rcout << "  -- split_min_events: " << split_min_events << std::endl;
    Rcout << "  -- pred_horizon: " << pred_horizon << std::endl;
    Rcout << std::endl << std::endl;
   }

   break;

  case TREE_CLASSIFICATION:

   forest = std::make_unique<ForestClassification>(data->n_cols_y);

   if(verbosity > 3){
    Rcout << "initializing classification forest" << std::endl;
    Rcout << "  -- n_class: " << data->n_cols_y << std::endl;
    Rcout << std::endl << std::endl;
   }

   break;

  case TREE_REGRESSION:

   forest = std::make_unique<ForestRegression>();

   if(verbosity > 3){
    Rcout << "initializing regression forest" << std::endl;
    Rcout << std::endl << std::endl;
   }

   break;

  default:

   Rcpp::stop("unrecognized tree type");
   break;

  }

  // does the forest need to be grown?
  bool grow_mode = loaded_forest.size() == 0;

  forest->init(std::move(data),
               tree_seeds,
               n_tree,
               mtry,
               sample_with_replacement,
               sample_fraction,
               grow_mode,
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
               pred_aggregate,
               pd_type,
               pd_x_vals,
               pd_x_cols,
               pd_probs,
               oobag,
               oobag_eval_type,
               oobag_eval_every,
               oobag_R_function,
               n_thread,
               verbosity);

   // Load forest object if it was already grown
   if(!grow_mode){

    uword n_obs = loaded_forest["n_obs"];

    std::vector<uvec>                rows_oobag   = loaded_forest["rows_oobag"];
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

     temp.load(n_tree, n_obs, rows_oobag, cutpoint, child_left,
               coef_values, coef_indices, leaf_pred_indx,
               leaf_pred_prob, leaf_pred_chaz, leaf_summary,
               pd_type, pd_x_vals, pd_x_cols, pd_probs);

    } else if (tree_type == TREE_CLASSIFICATION){

     std::vector<std::vector<vec>> leaf_pred_prob = loaded_forest["leaf_pred_prob"];

     auto& temp = dynamic_cast<ForestClassification&>(*forest);

     uword n_class = y.n_cols;

     temp.load(n_tree, n_obs, n_class, rows_oobag, cutpoint, child_left,
               coef_values, coef_indices, leaf_pred_prob, leaf_summary,
               pd_type, pd_x_vals, pd_x_cols, pd_probs);


    } else if (tree_type == TREE_REGRESSION){

     std::vector<std::vector<vec>> leaf_pred_prob = loaded_forest["leaf_pred_prob"];

     auto& temp = dynamic_cast<ForestRegression&>(*forest);

     temp.load(n_tree, n_obs, rows_oobag, cutpoint, child_left,
               coef_values, coef_indices, leaf_pred_prob, leaf_summary,
               pd_type, pd_x_vals, pd_x_cols, pd_probs);


    }

   }

   if(run_forest){ forest->run(oobag); }

   if(pred_mode){

    result.push_back(forest->get_predictions(), "pred_new");

   } else if (grow_mode) {

    if (oobag) result.push_back(forest->get_predictions(), "pred_oobag");

    List eval_oobag;
    eval_oobag.push_back(forest->get_oobag_eval(), "stat_values");
    eval_oobag.push_back(oobag_eval_type_R, "stat_type");
    result.push_back(eval_oobag, "eval_oobag");

   }

   if(write_forest){

    List forest_out;
    forest_out.push_back(n_obs, "n_obs");
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

    } else if (tree_type == TREE_CLASSIFICATION){

     auto& temp = dynamic_cast<ForestClassification&>(*forest);
     forest_out.push_back(temp.get_leaf_pred_prob(), "leaf_pred_prob");

    } else if (tree_type == TREE_REGRESSION){

     auto& temp = dynamic_cast<ForestRegression&>(*forest);
     forest_out.push_back(temp.get_leaf_pred_prob(), "leaf_pred_prob");

    }

    result.push_back(forest_out, "forest");

   }

   if(vi_type != VI_NONE){

    vec vi_output;
    if(run_forest){
     if(vi_type == VI_ANOVA){
      vi_output = forest->get_vi_numer() / forest->get_vi_denom();
     } else {
      vi_output = forest->get_vi_numer() / n_tree;
     }
    }
    result.push_back(vi_output, "importance");

   }

   if(pd_type != PD_NONE){
    result.push_back(forest->get_pd_values(), "pd_values");
   }

   return(result);

 }

