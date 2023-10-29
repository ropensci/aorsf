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

  TreeClassification(arma::uword n_class);

  TreeClassification(arma::uword n_obs,
                     arma::uvec& rows_oobag,
                     std::vector<double>& cutpoint,
                     std::vector<arma::uword>& child_left,
                     std::vector<arma::vec>& coef_values,
                     std::vector<arma::uvec>& coef_indices,
                     std::vector<arma::vec>& leaf_pred_prob,
                     std::vector<double>& leaf_summary);

  void resize_leaves(arma::uword new_size) override;

  double compute_split_score() override;

  void sprout_leaf_internal(arma::uword node_id) override;


  arma::uword n_class;

  // prob holds the predicted prob for each class
  std::vector<arma::vec> leaf_pred_prob;
  // summary (see Tree.h) holds class vote

 };

 } // namespace aorsf

#endif /* TREEClassification_H_ */
