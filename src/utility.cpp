/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "Utility.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 void print_mat(arma::mat& x,
                std::string label,
                arma::uword max_cols,
                arma::uword max_rows){

  uword nrow_print = max_rows-1;
  uword ncol_print = max_cols-1;

  if(x.n_rows < max_rows) nrow_print = x.n_rows-1;
  if(x.n_cols < max_cols) ncol_print = x.n_cols-1;


  Rcout << "   ---- view of " << label << " ---- " << std::endl << std::endl;
  Rcout << x.submat(0, 0, nrow_print, ncol_print);
  Rcout << std::endl << std::endl;

 }

 void print_vec(arma::vec& x,
                std::string label,
                arma::uword max_elem){

  uword n_print = max_elem-1;

  if(x.size() < n_print) n_print = x.size()-1;

  Rcout << "   ---- view of " << label << " ---- " << std::endl << std::endl;
  Rcout << x.subvec(0, n_print);
  Rcout << std::endl << std::endl;

 }

 void print_uvec(arma::uvec& x,
                 std::string label,
                 arma::uword max_elem){

  vec x_vec = conv_to<vec>::from(x);

  print_vec(x_vec, label, max_elem);

 }

 }
