/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef COXPH_H
#define COXPH_H

#include <armadillo>
#include "globals.h"


 namespace aorsf {

 // cholesky decomposition
 //
 // @description this function is copied from the survival package and
 //    translated into arma.
 //
 // @param vmat matrix with covariance estimates
 // @param n_vars the number of predictors used in the current node
 //
 // prepares vmat for cholesky_solve()


 void cholesky_decomp(arma::mat& vmat);

 // solve cholesky decomposition
 //
 // @description this function is copied from the survival package and
 //   translated into arma. Prepares u, the vector used to update beta.
 //
 // @param vmat matrix with covariance estimates
 // @param n_vars the number of predictors used in the current node
 //
 //
 void cholesky_solve(arma::mat& vmat,
                     arma::vec& u);

 // invert the cholesky in the lower triangle
 //
 // @description this function is copied from the survival package and
 //   translated into arma. Inverts vmat
 //
 // @param vmat matrix with covariance estimates
 // @param n_vars the number of predictors used in the current node
 //

 void cholesky_invert(arma::mat& vmat);

 // run the newton raphson procedure
 //
 // @description identify a linear combination of predictors.
 //   This function is copied from the survival package and
 //   translated into arma with light modifications for efficiency.
 //   The procedure works with the partial likelihood function
 //   of the Cox model. All inputs are described above
 //   in newtraph_cph_iter()
 //
 arma::mat coxph_fit(arma::mat& x_node,
                     arma::mat& y_node,
                     arma::vec& w_node,
                     arma::vec& XB,
                     bool do_scale,
                     int ties_method,
                     double epsilon,
                     arma::uword iter_max);

 }

#endif /* COXPH_H */

