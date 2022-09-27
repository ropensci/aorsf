/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef TREE_H_
#define TREE_H_

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

#include "Data.h"
#include "globals.h"

 namespace aorsf {

 class Tree {

 public:

  Tree() = default;

  Tree(Data* data,
       int leaf_min_obs,
       int mtry){

   this->data = data;

   arma::uword guess = std::ceil(
    0.5 * data->get_n_rows() / leaf_min_obs
   );

   coef.zeros(guess, mtry);
   coef_indices.zeros(guess, mtry);
   cutpoint.zeros(guess);
   next_left_node.zeros(guess);
   pred.zeros(guess, 1);
   pred_indices.zeros(guess, 3);


  };



  // @description sample weights to mimic a bootstrap sample
  // Note: the sampling extension for RcppArmadillo can only
  // be defined once. So, all functions that use sample need
  // to be defined in this file, unless we move the inclusion
  // of RcppArmadilloExtension/sample.h to another place.
  void draw_bootstrap_sample() {

   // s is the number of times you might get selected into
   // a bootstrap sample. Realistically this won't be >10,
   Rcpp::IntegerVector s = Rcpp::seq(0, 10);

   // compute probability of being selected into the bootstrap
   // 0 times, 1, times, ..., 9 times, or 10 times.

   arma::uword n_rows = data->get_n_rows();

   Rcpp::NumericVector probs = Rcpp::dbinom(s, n_rows, 1.0/n_rows, false);

   arma::vec boot_wts = Rcpp::as<arma::vec>(
    Rcpp::RcppArmadillo::sample(s, n_rows, true, probs)
   );

   if(data->has_weights()){

    boot_wts = boot_wts % data->w;

   }

   rows_inbag = arma::find(boot_wts);
   boot_wts = boot_wts(rows_inbag);

   if(VERBOSITY > 0){
    Rcpp::Rcout << "------------ in-bag rows ------------"   << std::endl;
    Rcpp::Rcout << rows_inbag                   << std::endl << std::endl;
    Rcpp::Rcout << "------------ boot weights -----------"   << std::endl;
    Rcpp::Rcout << boot_wts                     << std::endl << std::endl;

   }

  }

  void grow() {

   arma::mat x_inbag = data->x_rows(rows_inbag);
   arma::mat y_inbag = data->y_rows(rows_inbag);
   arma::vec w_inbag = data->w_subvec(rows_inbag);

   arma::vec node_assignments(rows_inbag.size(), arma::fill::zeros);
   arma::uvec nodes_to_grow(1, arma::fill::zeros);
   arma::uword nodes_max_true = 0;
   arma::uword leaf_node_counter = 0;
   arma::uword leaf_node_index_counter = 0;

   if(VERBOSITY > 0){

    arma::uword temp_uword_1, temp_uword_2;

    if(x_inbag.n_rows < 5)
     temp_uword_1 = x_inbag.n_rows-1;
    else
     temp_uword_1 = 5;

    if(x_inbag.n_cols < 5)
     temp_uword_2 = x_inbag.n_cols-1;
    else
     temp_uword_2 = 4;

    Rcpp::Rcout << "---- here is a view of x_inbag ---- " << std::endl;
    Rcpp::Rcout << x_inbag.submat(0, 0, temp_uword_1, temp_uword_2);
    Rcpp::Rcout << std::endl << std::endl;

   }

  };

  // INPUTS

  // Pointer to original data
  Data* data;

  // which rows of data are used to grow the tree
  arma::uvec rows_inbag;

  // coefficients for linear combinations;
  // one row per variable (mtry rows), one column per node
  // leaf nodes have all coefficients=0
  arma::mat coef;

  // indices of the predictors used by
  arma::umat coef_indices;

  // cutpoints used to split the node
  arma::vec cutpoint;

  // directions to the next node (right node = left node + 1)
  arma::uvec next_left_node;

  // predicted values (only in leaf nodes)
  arma::mat pred;

  // indices of predicted values for each leaf node
  arma::umat pred_indices;


 protected:

 };

 } // namespace aorsf

#endif /* TREE_H_ */
