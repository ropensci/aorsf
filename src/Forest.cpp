//  Forest.cpp

#include <RcppArmadillo.h>
#include "Forest.h"
#include "Tree.h"

using namespace arma;
using namespace Rcpp;

namespace aorsf {

Forest::Forest(){ }

void Forest::init(std::unique_ptr<Data> input_data,
                  Rcpp::IntegerVector& tree_seeds,
                  int n_tree,
                  int mtry,
                  VariableImportance vi_type,
                  double leaf_min_events,
                  double leaf_min_obs,
                  SplitRule split_rule,
                  double split_min_events,
                  double split_min_obs,
                  double split_min_stat,
                  int    split_max_retry,
                  LinearCombo lincomb_type,
                  double lincomb_eps,
                  int    lincomb_iter_max,
                  bool   lincomb_scale,
                  double lincomb_alpha,
                  int    lincomb_df_target,
                  PredType pred_type,
                  double pred_horizon,
                  bool   oobag_pred,
                  int    oobag_eval_every){

 this->data = std::move(input_data);
 this->tree_seeds = tree_seeds;
 this->n_tree = n_tree;
 this->mtry = mtry;
 this->vi_type = vi_type;
 this->leaf_min_events = leaf_min_events;
 this->leaf_min_obs = leaf_min_obs;
 this->split_rule = split_rule;
 this->split_min_events = split_min_events;
 this->split_min_obs = split_min_obs;
 this->split_min_stat = split_min_stat;
 this->split_max_retry = split_max_retry;
 this->lincomb_type = lincomb_type; this->lincomb_eps = lincomb_eps;
 this->lincomb_iter_max = lincomb_iter_max;
 this->lincomb_scale = lincomb_scale;
 this->lincomb_alpha = lincomb_alpha;
 this->lincomb_df_target = lincomb_df_target;
 this->pred_type = pred_type;
 this->pred_horizon = pred_horizon;
 this->oobag_pred = oobag_pred;
 this->oobag_eval_every = oobag_eval_every;

  if(VERBOSITY > 0){
  Rcout << "------------ dimensions ------------"   << std::endl;
  Rcout << "N obs total: "     << data->get_n_rows() << std::endl;
  Rcout << "N columns total: " << data->get_n_cols() << std::endl;
  Rcout << "------------------------------------";
  Rcout << std::endl << std::endl;
 }

 // sample weights to mimic a bootstrap sample
 this->bootstrap_select_times = seq(0, 10);

 uword n_rows = data->get_n_rows();

 // compute probability of being selected into the bootstrap
 // 0 times, 1, times, ..., 9 times, or 10 times.
 this->bootstrap_select_probs = dbinom(bootstrap_select_times,
                                       n_rows,
                                       1.0 / n_rows,
                                       false);

}

void Forest::grow(){ }

void Forest::run(){ }

}

