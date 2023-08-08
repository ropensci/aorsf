/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef TREE_H_
#define TREE_H_


#include "Data.h"
#include "globals.h"

 namespace aorsf {

 class Tree {

 public:

  Tree();

  // deleting the copy constructor
  Tree(const Tree&) = delete;
  // deleting the copy assignment operator
  Tree& operator=(const Tree&) = delete;

  void init(Data* data,
            double leaf_min_obs,
            double leaf_min_events,
            int mtry);

  void grow();

  // Pointer to original data
  const Data* data;

  // which node each inbag observation is currently in.
  arma::vec node_assignments;

  // which rows of data are held out while growing the tree
  arma::uvec rows_oobag;

  // coefficients for linear combinations;
  // one row per variable (mtry rows), one column per node
  // leaf nodes have all coefficients=0
  std::vector<arma::mat> coef_values;

  // indices of the predictors used by a node
  std::vector<arma::umat> coef_indices;

  // cutpoints used to split the nodes
  std::vector<double> cutpoint;

  // directions to the next node (right node = left node + 1)
  std::vector<arma::uword> child_left;

  // leaf values (only in leaf nodes)
  std::vector<arma::mat> leaf_values;

  // predicted values for out-of-bag rows
  arma::vec pred_oobag;

  // contains the node ID of each leaf node
  std::vector<size_t> leaf_index;

 protected:



 };

 } // namespace aorsf

#endif /* TREE_H_ */
