/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "Tree.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 Tree::Tree(){ }

 void Tree::init(Data* data,
                 double leaf_min_obs,
                 double leaf_min_events,
                 int mtry){


  arma::uword guess = std::ceil(
   0.5 * data->get_n_rows() / leaf_min_obs
  );

  // coef.zeros(guess, forest->mtry);
  // coef_indices.zeros(guess, forest->mtry);
  // cutpoint.zeros(guess);
  // next_left_node.zeros(guess);
  // leaf_values.zeros(guess, 1);
  // leaf_indices.zeros(guess, 3);

 }


 void Tree::grow(){


  // vec boot_wts = as<vec>(
  //  sample(forest->bootstrap_select_times,
  //         forest->data->get_n_rows(),
  //         true,
  //         forest->bootstrap_select_probs)
  // );

  //
  // if(forest->data->has_weights()) boot_wts = boot_wts % data->w;
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

