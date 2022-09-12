/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <armadillo>

#include "globals.h"
#include "Tree.h"

 namespace aorsf {

 Tree::Tree() :
   x_input(0),
   y_input(0),
   mtry(0),
   max_retry(DEFAULT_MAX_RETRY),
   split_rule(DEFAULT_SPLITRULE_SURVIVAL),
   n_split(DEFAULT_N_SPLIT),
   leaf_min_obs(0),
   split_min_obs(0),
   split_min_stat(DEFAULT_SPLIT_MIN_STAT),
   pred_type(DEFAULT_PRED_TYPE),
   oobag_eval_every(0),
   variable_importance(DEFAULT_IMPORTANCE),
   seed(0) {

   }

 void Tree::init(
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
   ) {


  this->x_input = x_input;
  this->y_input = y_input;
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

 }



 } // namespace aorsf
