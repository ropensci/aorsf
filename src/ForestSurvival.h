/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef FORESTSURVIVAL_H_
#define FORESTSURVIVAL_H_


#include "globals.h"
#include "Forest.h"
#include "TreeSurvival.h"

// [[Rcpp::depends(RcppArmadillo)]]

 namespace aorsf {

 class ForestSurvival: public Forest {
 public:

  // Default constructor with no arguments
  ForestSurvival();

  // Expect to redefine constructors for survival/classif/regression
  virtual ~ForestSurvival() = default;

  // Don't allow unitialized ForestSurvivals
  ForestSurvival(const ForestSurvival&) = delete;
  ForestSurvival& operator=(const ForestSurvival&) = delete;

 protected:

 };

 } // namespace aorsf

#endif /* FOREST_H_ */
