/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef DATA_H_
#define DATA_H_

#include <Rcpp.h>
#include <armadillo>
#include "globals.h"

 namespace aorsf {

 class Data {

 public:

  Data() = default;

  Data(arma::mat& x,
       arma::vec& y_dbl,
       arma::ivec& y_int,
       arma::vec& weights) {

   // arma::mat x_       = arma::mat(x.begin(), x.nrow(), x.ncol(), false);
   // arma::vec y_dbl_  = arma::vec(y_dbl.begin(), y_dbl.length(), false);
   // arma::ivec y_int_ = arma::ivec(y_int.begin(), y_int.length(), false);
   // arma::vec weights_ = arma::vec(weights.begin(), weights.length(), false);

   this->x = x;
   this->y_dbl = y_dbl;
   this->y_int = y_int;
   this->weights = weights;

  }

  arma::uword get_n_rows() const {
   return(x.n_rows);
  }

  arma::uword get_n_cols() const {
   return(x.n_cols);
  }

  arma::mat get_x_rows(arma::uvec vector_of_row_indices) const {
   return(x.rows(vector_of_row_indices));
  }

  arma::mat get_x_cols(arma::uvec vector_of_column_indices) const {
   return(x.cols(vector_of_column_indices));
  }

  double get_x_at(arma::uword i, arma::uword j) const {
   return(x.at(i, j));
  }

  arma::vec get_y_dbl_rows(arma::uvec vector_of_row_indices) const {
   return(y_dbl.rows(vector_of_row_indices));
  }

  arma::ivec get_y_int_rows(arma::uvec vector_of_row_indices) const {
   return(y_int.rows(vector_of_row_indices));
  }

  double get_y_dbl_at(arma::uword i) const {
   return(y_dbl.at(i));
  }

  int get_y_int_at(arma::uword i) const {
   return(y_int.at(i));
  }

 private:
  arma::mat x;
  arma::vec y_dbl;
  arma::ivec y_int;
  arma::vec weights;
 };


 } // namespace aorsf

#endif /* DATA_H_ */
