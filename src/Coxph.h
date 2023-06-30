/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef COXPH_H
#define COXPH_H

#include <armadillo>


 namespace aorsf {

 // scale observations in predictor matrix
 //
 // @description this scales inputs in the same way as
 //   the survival::coxph() function. The main reasons we do this
 //   are to avoid exponential overflow and to prevent the scale
 //   of inputs from impacting the estimated beta coefficients.
 //   E.g., you can try multiplying numeric inputs by 100 prior
 //   to calling orsf() with orsf_control_fast(do_scale = FALSE)
 //   and you will see that you get back a different forest.
 //
 // @param x_node matrix of predictors
 // @param w_node replication weights
 // @param x_transforms matrix used to store the means and scales
 //
 // @return modified x_node and x_transform filled with values
 //
 arma::mat coxph_scale(arma::mat& x_node,
                       arma::vec& w_node);

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
 arma::vec newtraph_cph(arma::mat& x_node,
                        arma::mat& y_node,
                        arma::vec& w_node,
                        int ties_method,
                        double cph_eps,
                        arma::uword cph_iter_max,
                        char oobag_importance_type);

}

#endif /* COXPH_H */

