/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <armadillo>
#include <RcppArmadillo.h>
 // [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include "globals.h"
#include "Tree.h"

 namespace aorsf {

 void Tree::guess_max_nodes() {

  // arma::uword guess = std::ceil(
  //  0.5 * data->get_n_rows() / this->leaf_min_obs
  // );

  // int a = 2;
  // int b = 4;
  //
  // this->coef = Rcpp::NumericMatrix(a,b);
  // this->coef_indices = arma::umat(a, b);

  // cutpoint.zeros(guess);
  // left_node.zeros(guess);
  // pred.zeros(guess, 1);

 }

 Tree::Tree() :
   data(0),
   mtry(0),
   max_retry(DEFAULT_MAX_RETRY),
   split_rule(DEFAULT_SPLITRULE_SURVIVAL),
   n_split(DEFAULT_N_SPLIT),
   leaf_min_obs(0),
   split_min_obs(0),
   split_min_stat(DEFAULT_SPLIT_MIN_STAT),
   pred_type(DEFAULT_PRED_TYPE),
   oobag_eval_every(0),
   variable_importance(DEFAULT_IMPORTANCE),
   seed(0) {

   }

 } // namespace aorsf
