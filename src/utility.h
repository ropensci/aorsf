/*-----------------------------------------------------------------------------
 This file is part of aorsf.
Author: Byron C Jaeger
aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef UTILITY_H
#define UTILITY_H

#include <armadillo>
#include <Rcpp.h>


 namespace aorsf {

 void print_mat(arma::mat& x,
                std::string label,
                arma::uword max_cols,
                arma::uword max_rows);

 void print_vec(arma::vec& x,
                std::string label,
                arma::uword max_elem);

 void print_uvec(arma::uvec& x,
                 std::string label,
                 arma::uword max_elem);

 }

#endif /* UTILITY_H */