/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef DATA_H_
#define DATA_H_

#include <armadillo>
#include "globals.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 class Data {

 public:

  Data() = default;

  Data(arma::mat& x,
       arma::mat& y,
       arma::vec& w) {

   this->x = x;
   this->y = y;
   this->w = w;

   this->n_rows = x.n_rows;
   this->n_cols = x.n_cols;
   this->has_weights = !w.empty();

  }

  Data(const Data&) = delete;
  Data& operator=(const Data&) = delete;

  arma::uword get_n_rows() {
   return(n_rows);
  }

  arma::uword get_n_cols() {
   return(n_cols);
  }

  arma::mat x_rows(arma::uvec& vector_of_row_indices) {
   return(x.rows(vector_of_row_indices));
  }

  arma::mat x_cols(arma::uvec& vector_of_column_indices) {
   return(x.cols(vector_of_column_indices));
  }

  arma::mat y_rows(arma::uvec& vector_of_row_indices) {
   return(y.rows(vector_of_row_indices));
  }

  arma::mat y_cols(arma::uvec& vector_of_column_indices) {
   return(y.cols(vector_of_column_indices));
  }

  arma::mat x_submat(arma::uvec& vector_of_row_indices,
                     arma::uvec& vector_of_column_indices){
   return(x.submat(vector_of_row_indices, vector_of_column_indices));
  }

  arma::mat y_submat(arma::uvec& vector_of_row_indices,
                     arma::uvec& vector_of_column_indices){
   return(y.submat(vector_of_row_indices, vector_of_column_indices));
  }

  arma::vec w_subvec(arma::uvec& vector_of_indices){
   return(w(vector_of_indices));
  }

  // member variables

  arma::uword n_cols;
  arma::uword n_rows;
  arma::vec w;

  bool has_weights;

 private:

  arma::mat x;
  arma::mat y;

 };


 } // namespace aorsf

#endif /* DATA_H_ */
