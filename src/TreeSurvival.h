/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef TREESURVIVAL_H_
#define TREESURVIVAL_H_

#include "globals.h"
#include "Tree.h"

// [[Rcpp::depends(RcppArmadillo)]]

 namespace aorsf {

 class TreeSurvival : public Tree {
 public:

  TreeSurvival();

  // Don't allow unitialized trees
  TreeSurvival(const TreeSurvival&) = delete;
  TreeSurvival& operator=(const TreeSurvival&) = delete;

  virtual ~TreeSurvival() override = default;

 };

 } // namespace aorsf

#endif /* TREESURVIVAL_H_ */
