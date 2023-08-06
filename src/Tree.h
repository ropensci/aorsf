/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef TREE_H_
#define TREE_H_


#include "Data.h"
#include "globals.h"
#include "Forest.h"

 namespace aorsf {

 class Tree {

 public:

  // Construct trees from an existing forest
  Tree(Forest* forest);

  // deleting the copy constructor
  Tree(const Tree&) = delete;
  // deleting the copy assignment operator
  Tree& operator=(const Tree&) = delete;

  void grow();

  // which rows of data are held out while growing the tree
  arma::uvec rows_oobag;

  // coefficients for linear combinations;
  // one row per variable (mtry rows), one column per node
  // leaf nodes have all coefficients=0
  arma::mat coef;

  // indices of the predictors used by a node
  arma::umat coef_indices;

  // cutpoints used to split the nodes
  arma::vec cutpoint;

  // directions to the next node (right node = left node + 1)
  arma::uvec next_left_node;

  // leaf values (only in leaf nodes)
  arma::mat leaf_values;

  // predicted values for out-of-bag rows
  arma::vec pred_oobag;

  // indices of predicted values for each leaf node
  arma::umat leaf_indices;

 protected:

  // Pointer to a forest
  Forest* forest;


 };

 } // namespace aorsf

#endif /* TREE_H_ */
