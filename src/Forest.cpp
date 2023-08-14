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
                  arma::uword n_tree,
                  arma::uword mtry,
                  VariableImportance vi_type,
                  double leaf_min_events,
                  double leaf_min_obs,
                  SplitRule split_rule,
                  double split_min_events,
                  double split_min_obs,
                  double split_min_stat,
                  arma::uword split_max_cuts,
                  arma::uword split_max_retry,
                  LinearCombo lincomb_type,
                  double lincomb_eps,
                  arma::uword lincomb_iter_max,
                  bool   lincomb_scale,
                  double lincomb_alpha,
                  arma::uword lincomb_df_target,
                  PredType pred_type,
                  bool   pred_mode,
                  double pred_horizon,
                  bool   oobag_pred,
                  arma::uword oobag_eval_every){

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
 this->split_max_cuts = split_max_cuts;
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
  Rcout << "------------ input data dimensions ------------"   << std::endl;
  Rcout << "N obs total: "     << data->get_n_rows() << std::endl;
  Rcout << "N columns total: " << data->get_n_cols() << std::endl;
  Rcout << "-----------------------------------------------";
  Rcout << std::endl << std::endl;
 }

}

// growInternal() in ranger
void Forest::plant() {

 trees.reserve(n_tree);

 for (arma::uword i = 0; i < n_tree; ++i) {
  trees.push_back(std::make_unique<Tree>());
 }

}

void Forest::grow(Function& lincomb_R_function){


 for(uword i = 0; i < n_tree; ++i){

  trees[i]->init(data.get(),
                 tree_seeds[i],
                 mtry,
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
                 lincomb_df_target);

  trees[i]->grow();


 }

 double x_dbl = 1.0;

 NumericMatrix test_mat = lincomb_R_function(x_dbl);

 arma::mat test_mat_arma(test_mat.begin(),
                         test_mat.nrow(),
                         test_mat.ncol(), false);

 Rcout << "--- test R function output ---" << std::endl << std::endl;
 Rcout << test_mat_arma << std::endl;

 // result.push_back(test_mat_arma, "test");



}

void Forest::run(){ }

}

