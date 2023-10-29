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

 namespace aorsf {

 class Data {

 public:

  Data() = default;

  Data(const Data&) = delete;
  Data& operator=(const Data&) = delete;

  virtual ~Data() = default;

  Data(arma::mat& x,
       arma::mat& y,
       arma::vec& w) {

   this->x = x;
   this->y = y;
   this->w = w;

   this->n_rows = x.n_rows;
   this->n_cols_x = x.n_cols;
   this->n_cols_y = y.n_cols;
   this->has_weights = !w.empty();
   this->saved_values.resize(n_cols_x);

  }


  arma::uword get_n_rows() {
   return(n_rows);
  }

  arma::uword get_n_cols_x() {
   return(n_cols_x);
  }

  arma::mat& get_x(){
   return(x);
  }

  arma::mat& get_y(){
   return(y);
  }

  arma::vec& get_w(){
   return(w);
  }

  arma::mat x_rows(arma::uvec& row_indices) {
   return(x.rows(row_indices));
  }

  arma::mat x_cols(arma::uvec& column_indices) {
   return(x.cols(column_indices));
  }

  arma::mat y_rows(arma::uvec& row_indices) {
   return(y.rows(row_indices));
  }

  arma::mat y_cols(arma::uvec& column_indices) {
   return(y.cols(column_indices));
  }

  arma::mat x_submat(arma::uvec& row_indices,
                     arma::uvec& column_indices){
   return(x.submat(row_indices, column_indices));
  }

  arma::mat y_submat(arma::uvec& row_indices,
                     arma::uvec& column_indices){
   return(y.submat(row_indices, column_indices));
  }

  arma::vec w_subvec(arma::uvec& indices){
   return(w(indices));
  }

  // multiply X matrix by lincomb coefficients
  // without taking a sub-matrix of X
  arma::vec x_submat_mult_beta(arma::uvec& x_rows,
                               arma::uvec& x_cols,
                               arma::vec& beta){

   arma::vec out(x_rows.size());
   arma::uword i = 0;

   for(auto row : x_rows){
    arma::uword j = 0;
    for(auto col : x_cols){
     out[i] += x.at(row, col) * beta[j];
     j++;
    }
    i++;
   }

   return(out);

  }

  void permute_col(arma::uword j, std::mt19937_64& rng){

   arma::vec x_j = x.unsafe_col(j);
   // make and store a copy
   this->saved_values[j] = arma::vec(x_j.begin(), x_j.size(), true);
   // shuffle the vector in-place
   std::shuffle(x_j.begin(), x_j.end(), rng);

  }

  void save_col(arma::uword j){
   saved_values[j] = x.col(j);
  }

  void restore_col(arma::uword j){
   x.col(j) = saved_values[j];
  }

  void fill_col(double value, uword j){
    x.col(j).fill(value);
  }


  // member variables

  arma::uword n_cols_x;
  arma::uword n_cols_y;
  arma::uword n_rows;
  arma::vec w;

  // for multi-column ops (e.g., partial dependence)
  std::vector<arma::vec> saved_values;

  bool has_weights;

 private:

  arma::mat x;
  arma::mat y;

 };


 } // namespace aorsf

#endif /* DATA_H_ */
