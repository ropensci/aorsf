/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <armadillo>

 namespace aorsf {

 typedef unsigned int uint;

 template <typename T>
 using svec = std::vector<T>;

 // Tree types
 enum TreeType {
  TREE_CLASSIFICATION = 1,
  TREE_REGRESSION = 2,
  TREE_SURVIVAL = 3,
  TREE_PROBABILITY = 4
 };

 // Variable importance
 enum VariableImportance {
  VI_NONE = 0,
  VI_NEGATE = 1,
  VI_PERMUTE = 2,
  VI_ANOVA = 3
 };

 // Split mode
 enum SplitRule {
  SPLIT_LOGRANK = 1,
  SPLIT_CONCORD = 2
 };

 // Linear combination method
 enum LinearCombo {
  NEWTON_RAPHSON = 1,
  RANDOM_COEFS = 2,
  R_FUNCTION = 3
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
 const int DEFAULT_N_THREADS = 1;

 const VariableImportance DEFAULT_IMPORTANCE = VI_NONE;

 const double DEFAULT_SPLIT_MAX_RETRY = 1;


 const double DEFAULT_LEAF_MIN_EVENTS = 1;
 const double DEFAULT_LEAF_MIN_OBS = 5;

 const SplitRule DEFAULT_SPLITRULE = SPLIT_LOGRANK;
 const double    DEFAULT_SPLIT_MIN_EVENTS = 5;
 const double    DEFAULT_SPLIT_MIN_OBS = 10;
 const double    DEFAULT_SPLIT_MIN_STAT = 3.84;

 const arma::uword DEFAULT_SPLIT_MAX_CUTS = 5;
 const arma::uword DEFAULT_MAX_RETRY = 3;

 const LinearCombo DEFAULT_LINCOMB = NEWTON_RAPHSON;
 const double      DEFAULT_LINCOMB_EPS = 1e-9;
 const arma::uword DEFAULT_LINCOMB_ITER_MAX = 20;
 const bool        DEFAULT_LINCOMB_SCALE = true;
 const double      DEFAULT_LINCOMB_ALPHA = 0.5;
 const arma::uword DEFAULT_LINCOMB_TIES_METHOD = 1;

 const double DEFAULT_ANOVA_VI_PVALUE = 0.01;

 const PredType DEFAULT_PRED_TYPE = RISK;

 const int VERBOSITY = 0;

 // Interval to print progress in seconds
 const double STATUS_INTERVAL = 15.0;


 } // namespace aorsf

#endif /* GLOBALS_H_ */
