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

  TreeSurvival(double leaf_min_events,
               double split_min_events,
               arma::vec* unique_event_times,
               arma::vec* pred_horizon);

  TreeSurvival(std::vector<double>& cutpoint,
               std::vector<arma::uword>& child_left,
               std::vector<arma::vec>& coef_values,
               std::vector<arma::uvec>& coef_indices,
               std::vector<arma::vec>& leaf_pred_indx,
               std::vector<arma::vec>& leaf_pred_prob,
               std::vector<arma::vec>& leaf_pred_chaz,
               std::vector<double>& leaf_summary,
               arma::vec* pred_horizon);

  double compute_max_leaves() override;

  void resize_leaves(arma::uword new_size) override;

  bool is_col_splittable(arma::uword j) override;

  bool is_node_splittable_internal() override;

  arma::uvec find_cutpoints() override;

  double compute_split_score() override;

  double score_logrank();

  double compute_mortality(arma::mat& leaf_data);

  void sprout_leaf(uword node_id) override;

  void predict_value(arma::mat* pred_output,
                     arma::vec* pred_denom,
                     PredType pred_type,
                     bool oobag) override;

  std::vector<arma::vec>& get_leaf_pred_indx(){
   return(leaf_pred_indx);
  }

  std::vector<arma::vec>& get_leaf_pred_prob(){
   return(leaf_pred_prob);
  }

  std::vector<arma::vec>& get_leaf_pred_chaz(){
   return(leaf_pred_chaz);
  }

  double compute_prediction_accuracy(arma::vec& preds) override;

  std::vector<arma::vec> leaf_pred_indx;
  std::vector<arma::vec> leaf_pred_prob;
  std::vector<arma::vec> leaf_pred_chaz;

  // pointer to event times in forest
  arma::vec* unique_event_times;

  // prediction times
  arma::vec* pred_horizon;

  double leaf_min_events;
  double split_min_events;

 };

 } // namespace aorsf

#endif /* TREESURVIVAL_H_ */
