/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef NODESPLITSTATS_H
#define NODESPLITSTATS_H

#include <armadillo>
#include <Rcpp.h>


 namespace aorsf {

 arma::uvec node_find_cps(const arma::mat& y_node,
                          const arma::vec& w_node,
                          const arma::vec& XB,
                          arma::uvec& XB_sorted,
                          double leaf_min_events,
                          double leaf_min_obs);

 void node_fill_group(arma::vec& group,
                      const arma::uvec& XB_sorted,
                      const arma::uword start,
                      const arma::uword stop,
                      const double value);


 double node_compute_lrt(arma::mat& y_node,
                         arma::vec& w_node,
                         arma::vec& group);

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
 Rcpp::List lrt_multi(arma::mat& y_node,
                      arma::mat& w_node,
                      arma::vec& XB,
                      arma::uword n_split,
                      double split_min_stat,
                      double leaf_min_events,
                      double leaf_min_obs);


 }

#endif /* NODESPLITSTATS_H */

