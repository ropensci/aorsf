/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef FOREST_H_
#define FOREST_H_


#include "Data.h"
#include "globals.h"

// [[Rcpp::depends(RcppArmadillo)]]

 namespace aorsf {

 class Forest {
 public:

  // Default constructor with no arguments
  Forest();

  // Expect to redefine constructors for survival/classif/regression
  virtual ~Forest() = default;

  // Don't allow unitialized Forests
  Forest(const Forest&) = delete;
  Forest& operator=(const Forest&) = delete;

  void init(
    const Data* data,
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
  ) {

   this->data = data;
   this->mtry = mtry;
   this->max_retry = max_retry;
   this->split_rule = split_rule;
   this->n_split = n_split;
   this->leaf_min_obs = leaf_min_obs;
   this->split_min_obs = split_min_obs;
   this->split_min_stat = split_min_stat;
   this->pred_type = pred_type;
   this->oobag_eval_every = oobag_eval_every;
   this->variable_importance = variable_importance;

  };

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
  const Data* data;

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

 };

 } // namespace aorsf

#endif /* FOREST_H_ */
