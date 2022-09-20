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



 } // namespace aorsf
