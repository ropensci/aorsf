//  Forest.cpp

#include <RcppArmadillo.h>
#include "Forest.h"
#include "Tree.h"

using namespace arma;
using namespace Rcpp;

namespace aorsf {

Forest::Forest(){ }

void Forest::init(std::unique_ptr<Data> input_data,
                  int n_tree,
                  Rcpp::IntegerVector& tree_seeds,
                  Rcpp::List& tree_params){

 this->data = std::move(input_data);

 uword n_rows = data->get_n_rows();

 // this->n_tree = n_tree;
 // this->tree_seeds = tree_seeds;
 // this->tree_objects = Rcpp::List(n_tree);

 // this->n_split = tree_params["n_split"];
 // this->mtry = tree_params["mtry"];
 // this->leaf_min_events = tree_params["leaf_min_events"];
 // this->leaf_min_obs = tree_params["leaf_min_obs"];

 // this->split_min_events = params["split_min_events"];
 // this->split_min_obs = params["split_min_obs"];
 // this->split_min_stat = params["split_min_stat"];
 // this->cph_method = params["cph_method"];
 // this->cph_eps = params["cph_eps"];
 // this->cph_iter_max = params["cph_iter_max"];
 // this->cph_do_scale = params["cph_do_scale"];
 // this->net_alpha = params["net_alpha"];
 // this->net_df_target = params["net_df_target"];
 // this->oobag_pred = params["oobag_pred"];
 // this->oobag_pred_type = params["oobag_pred_type"];
 // this->oobag_pred_horizon = params["oobag_pred_horizon"];
 // this->oobag_eval_every = params["oobag_eval_every"];
 // this->oobag_importance = params["oobag_importance"];
 // this->oobag_importance_type = params["oobag_importance_type"];
 // this->max_retry = params["max_retry"];
 // this->type_beta = params["type_beta"];


 if(VERBOSITY > 0){
  Rcout << "------------ dimensions ------------"   << std::endl;
  Rcout << "N obs total: "     << data->get_n_rows() << std::endl;
  Rcout << "N columns total: " << data->get_n_cols() << std::endl;
  Rcout << "------------------------------------";
  Rcout << std::endl << std::endl;
 }
 //
 // // sample weights to mimic a bootstrap sample
 //
 // // s is the number of times you might get selected into
 // // a bootstrap sample. Realistically this won't be >10,
 // // but it could technically be as big as n_row.
 this->bootstrap_select_times = seq(0, 10);

 // compute probability of being selected into the bootstrap
 // 0 times, 1, times, ..., 9 times, or 10 times.
 this->bootstrap_select_probs = dbinom(bootstrap_select_times,
                                       n_rows,
                                       1.0 / n_rows,
                                       false);

}

}

