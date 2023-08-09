/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "Tree.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 Tree::Tree(){ }

 void Tree::init(Data*  data,
                 int    seed,
                 int    mtry,
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
                 Rcpp::IntegerVector* bootstrap_select_times,
                 Rcpp::NumericVector* bootstrap_select_probs){


  this->data = data;
  this->seed = seed;
  this->mtry = mtry;
  this->leaf_min_events = leaf_min_events;
  this->leaf_min_obs = leaf_min_obs;
  this->split_rule = split_rule;
  this->split_min_events = split_min_events;
  this->split_min_obs = split_min_obs;
  this->split_min_stat = split_min_stat;
  this->split_max_retry = split_max_retry;
  this->lincomb_type = lincomb_type;
  this->lincomb_eps = lincomb_eps;
  this->lincomb_iter_max = lincomb_iter_max;
  this->lincomb_scale = lincomb_scale;
  this->lincomb_alpha = lincomb_alpha;
  this->lincomb_df_target = lincomb_df_target;

  // node_assignments.zeros(x_inbag.n_rows);
  nodes_to_grow.zeros(1);

  // leaf_node_counter = 0;
  // leaf_node_index_counter = 0;

 }

 void Tree::bootstrap(){

  // Initialize random number generator and set seed
  random_number_generator.seed(seed);

  uword i, draw, n = data->n_rows;

  // Start with all samples OOB
  vec boot_wts(n, fill::zeros);

  std::uniform_int_distribution<uword> unif_dist(0, n - 1);

  // Draw num_samples samples with replacement (num_samples_inbag out of n) as inbag and mark as not OOB
  for (i = 0; i < n; ++i) {
   draw = unif_dist(random_number_generator);
   ++boot_wts[draw];
  }

  // multiply boot_wts by user specified weights.
  if(data->has_weights){
   boot_wts = boot_wts = boot_wts % data->w;
  }

  uvec rows_inbag = find(boot_wts > 0);

  this->rows_oobag = find(boot_wts == 0);
  this->x_inbag = data->x_rows(rows_inbag);
  this->y_inbag = data->y_rows(rows_inbag);
  this->w_inbag = data->w_subvec(rows_inbag);

 }

 void Tree::grow(){


  bootstrap();

  Rcout << x_inbag << std::endl;

  //
  // arma::uvec rows_inbag = arma::find(boot_wts);
  //
  // // initalize the oob rows if oob predictions are being computed
  // if(oobag_pred)
  //  rows_oobag = arma::find(boot_wts != 0);
  //
  // boot_wts = boot_wts(rows_inbag);
  //
  // if(VERBOSITY > 0){
  //  Rcpp::Rcout << "------------ in-bag rows ------------"   << std::endl;
  //  Rcpp::Rcout << rows_inbag                   << std::endl << std::endl;
  //  Rcpp::Rcout << "------------ boot weights -----------"   << std::endl;
  //  Rcpp::Rcout << boot_wts                     << std::endl << std::endl;
  //
  // }
  //
  // arma::mat x_inbag = data->x_rows(rows_inbag);
  // arma::mat y_inbag = data->y_rows(rows_inbag);
  // arma::vec w_inbag = data->w_subvec(rows_inbag);
  //
  //
  //
  // // once the sub-matrix views are created, we do not use inbag rows
  // // (if we really need them, we can get them from oobag rows)
  // rows_inbag.resize(0);
  //
  // arma::vec node_assignments(rows_inbag.size(), arma::fill::zeros);
  // arma::uvec nodes_to_grow(1, arma::fill::zeros);
  // arma::uword nodes_max_true = 0;
  // arma::uword leaf_node_counter = 0;
  // arma::uword leaf_node_index_counter = 0;
  //
  //
  //
  // if(VERBOSITY > 0){
  //
  //  arma::uword temp_uword_1, temp_uword_2;
  //
  //  if(x_inbag.n_rows < 5)
  //   temp_uword_1 = x_inbag.n_rows-1;
  //  else
  //   temp_uword_1 = 5;
  //
  //  if(x_inbag.n_cols < 5)
  //   temp_uword_2 = x_inbag.n_cols-1;
  //  else
  //   temp_uword_2 = 4;
  //
  //  Rcpp::Rcout << "---- here is a view of x_inbag ---- " << std::endl;
  //  Rcpp::Rcout << x_inbag.submat(0, 0, temp_uword_1, temp_uword_2);
  //  Rcpp::Rcout << std::endl << std::endl;
  //
  // }



 } // Tree::grow

 } // namespace aorsf

