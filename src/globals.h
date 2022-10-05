/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <armadillo>

 namespace aorsf {

 // Tree types
 enum TreeType {
  TREE_CLASSIFICATION = 1,
  TREE_REGRESSION = 3,
  TREE_SURVIVAL = 5,
  TREE_PROBABILITY = 9
 };

 // Variable importance
 enum VariableImportance {
  VI_NONE = 0,
  VI_NEGATE = 1,
  VI_PERM_BREIMAN = 2,
  VI_ANOVA = 3
 };

 // Split mode
 enum SplitRule {
  SPLIT_LOGRANK = 1,
  SPLIT_GINI = 2
 };

 // Linear combination method
 enum LinearCombination {
  NEWTON_RAPHSON = 1,
  USER_FUNCTION = 2
 };

 // Prediction type
 enum PredType {
  RISK = 1,
  SURVIVAL = 2,
  CUMULATIVE_HAZARD = 3,
  MORTALITY = 4,
  MEAN = 5,
  PROBABILITY = 6,
  CLASS = 7,
  TERMINAL_NODES = 8
 };

 // Default values
 const int DEFAULT_N_TREE = 500;
 const int DEFAULT_N_THREADS = 0;
 const VariableImportance DEFAULT_IMPORTANCE = VI_NONE;

 const int DEFAULT_LEAF_MIN_OBS_CLASSIFICATION = 1;
 const int DEFAULT_LEAF_MIN_OBS_REGRESSION = 5;
 const int DEFAULT_LEAF_MIN_OBS_SURVIVAL = 10;
 const int DEFAULT_LEAF_MIN_OBS_PROBABILITY = 10;

 const int DEFAULT_MAX_RETRY = 3;

 const double DEFAULT_SPLIT_MIN_STAT = 3.84;
 const int DEFAULT_SPLIT_MIN_OBS = 10;

 const SplitRule DEFAULT_SPLITRULE_SURVIVAL = SPLIT_LOGRANK;

 const PredType DEFAULT_PRED_TYPE = RISK;
 const int DEFAULT_N_SPLIT = 5;

 const int VERBOSITY = 1;

 } // namespace aorsf

#endif /* GLOBALS_H_ */
