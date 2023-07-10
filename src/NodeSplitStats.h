/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef NODESPLITSTATS_H
#define NODESPLITSTATS_H

#include <armadillo>


 namespace aorsf {

 // Log rank test w/multiple cutpoints
 //
 // this function returns a cutpoint obtaining a local maximum
 // of the log-rank test (lrt) statistic. The default value (+Inf)
 // is really for diagnostic purposes. Put another way, if the
 // return value is +Inf (an impossible value for a cutpoint),
 // that means that we didn't find any valid cut-points and
 // the node cannot be grown with the current XB.
 //
 // if there is a valid cut-point, then the main side effect
 // of this function is to modify the group vector, which
 // will be used to assign observations to the two new nodes.
 //
 // @param group the vector that determines which node to send each
 //   observation to (left node = 0, right node = 1)
 // @param y_node matrix of outcomes
 // @param w_node vector of weights
 // @param XB linear combination of predictors
 //
 // the group vector is modified by this function and the value returned
 // is the maximal log-rank statistic across all the possible cutpoints.
 double lrt_multi(arma::mat& x_node,
                  arma::mat& y_node,
                  arma::mat& w_node,
                  arma::vec& XB,
                  arma::uword n_split,
                  double split_min_stat,
                  double leaf_min_obs,
                  double leaf_min_events);


 }

#endif /* NODESPLITSTATS_H */

