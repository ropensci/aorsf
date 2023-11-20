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
                     arma::uword n_class,
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

  arma::uword predict_value_internal(arma::uvec& pred_leaf_sort,
                                     arma::mat& pred_output,
                                     arma::vec& pred_denom,
                                     PredType pred_type,
                                     bool oobag) override;

  arma::uword find_safe_mtry() override;
  arma::uword find_safe_mtry_binary();
  arma::uword find_safe_mtry_multiclass();

  double compute_prediction_accuracy_internal(arma::mat& preds) override;

  arma::mat glm_fit() override;
  arma::mat glmnet_fit() override;
  arma::mat user_fit() override;

  uword get_n_col_vi() override;

  void fill_pred_values_vi(arma::mat& pred_values) override;

  std::vector<arma::vec>& get_leaf_pred_prob(){
   return(leaf_pred_prob);
  }

  arma::uword n_class;
  bool binary;

  arma::uvec splittable_y_cols;
  arma::uword y_col_split;

  // prob holds the predicted prob for each class
  std::vector<arma::vec> leaf_pred_prob;
  // summary (see Tree.h) holds class vote

 };

 } // namespace aorsf

#endif /* TREEClassification_H_ */
