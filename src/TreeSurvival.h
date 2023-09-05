/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef TREESURVIVAL_H_
#define TREESURVIVAL_H_


#include "Data.h"
#include "globals.h"
#include "Tree.h"

 namespace aorsf {

 class TreeSurvival: public Tree {

 public:

  TreeSurvival();
  TreeSurvival(const TreeSurvival&) = delete;
  TreeSurvival& operator=(const TreeSurvival&) = delete;

  double compute_max_leaves() override;

  bool is_node_splittable_internal() override;

  arma::uvec find_cutpoints() override;

 };

 } // namespace aorsf

#endif /* TREESURVIVAL_H_ */
