/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef TREE_H_
#define TREE_H_

#include <armadillo>

#include "globals.h"

 namespace aorsf {

 class Tree {
 public:

  // Default constructor with no arguments
  Tree();

  // Create from loaded forest
  // Tree(arma::mat& coef,
  //      arma::umat& coef_indices,
  //      arma::vec& cutpoint,
  //      arma::uvec& next_left_node,
  //      arma::mat& pred,
  //      arma::umat& pred_indices); // TODO: pred_indices to survival tree

  // Expect to redefine constructors for survival/classif/regression
  virtual ~Tree() = default;

  // Don't allow unitialized trees
  // Tree(const Tree&) = delete;
  // Tree& operator=(const Tree&) = delete;

  void init(
    double* x_input,
    double* y_input,
    const int mtry,
    const int max_retry,
    SplitRule split_rule,
    const int n_split,
    const int leaf_min_obs,
    const int split_min_obs,
    const int split_min_stat,
    PredType pred_type,
    const int oobag_eval_every,
    VariableImportance variable_importance,
    const int seed
  );

  int get_mtry() const {
   return mtry;
  }

  int get_max_retry() const {
   return max_retry;
  }

  SplitRule get_split_rule() const {
   return split_rule;
  }

  int get_n_split() const {
   return n_split;
  }

  int get_leaf_min_obs() const {
   return leaf_min_obs;
  }

  int get_split_min_obs() const {
   return split_min_obs;
  }

  int get_split_min_stat() const {
   return split_min_stat;
  }

  int get_pred_type() const {
   return pred_type;
  }

  int get_oobag_eval_every() const {
   return oobag_eval_every;
  }

  VariableImportance get_variable_importance() const {
   return variable_importance;
  }

  // INPUTS

  // Pointer to original data
  double* x_input;
  double* y_input;

  // number of predictors used to split a node
  int mtry;

  // maximum number of retry attempts to split a node
  int max_retry;

  // how to measure quality of a node split
  SplitRule split_rule;

  // number of cutpoints to assess during node split
  int n_split;

  // minimum number of observations needed in a leaf node
  int leaf_min_obs;

  // minimum number of observations needed to split a node
  int split_min_obs;

  // minimum value of split statistic needed to split a node
  int split_min_stat;

  // what type of oobag prediction to compute
  PredType pred_type;

  // evaluate oobag error every X trees
  int oobag_eval_every;

  // what type of variable importance to compute
  VariableImportance variable_importance;

  // random seed to be set before growing
  int seed;

 protected:

  // OUTPUTS

  // coefficients for linear combinations;
  // one row per variable (mtry rows), one column per node
  // leaf nodes have all coefficients=0
  arma::mat coef;

  // indices of the predictors used by
  arma::umat coef_indices;

  // cutpoints used to split the node
  arma::vec cutpoint;

  // directions to the next node (right node = left node + 1)
  arma::uvec next_left_node;

  // predicted values (only in leaf nodes)
  arma::mat pred;

  // indices of predicted values for each leaf node
  arma::umat pred_indices;

 };

 } // namespace aorsf

#endif /* TREE_H_ */
