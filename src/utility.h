/*-----------------------------------------------------------------------------
 This file is part of aorsf.
Author: Byron C Jaeger
aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef UTILITY_H
#define UTILITY_H

#include <armadillo>
#include <Rcpp.h>
#include "globals.h"


 namespace aorsf {

 /**
  * Split sequence start..end in num_parts parts with sizes as equal as possible.
  * @param result Result vector of size num_parts+1. Ranges for the parts are then result[0]..result[1]-1, result[1]..result[2]-1, ..
  * @param start minimum value
  * @param end maximum value
  * @param num_parts number of parts
  * @note: this function is taken directly from ranger.
  */
 void equalSplit(std::vector<uint>& result, uint start, uint end, uint num_parts);


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
