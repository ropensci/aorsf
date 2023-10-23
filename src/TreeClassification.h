/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef TREECLASSIFICATION_H_
#define TREECLASSIFICATION_H_


#include "Data.h"
#include "globals.h"
#include "Tree.h"

 namespace aorsf {

 class TreeClassification: public Tree {

 public:

  TreeClassification();

  TreeClassification(const TreeClassification&) = delete;
  TreeClassification& operator=(const TreeClassification&) = delete;

  virtual ~TreeClassification() override = default;

 };

 } // namespace aorsf

#endif /* TREEClassification_H_ */
