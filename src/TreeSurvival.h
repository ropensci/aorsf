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

  virtual ~TreeSurvival() override = default;

  TreeSurvival(double leaf_min_events,
               double split_min_events,
               arma::vec* unique_event_times,
               arma::vec* pred_horizon);

  TreeSurvival(arma::uword n_obs,
               arma::uvec& rows_oobag,
               std::vector<double>& cutpoint,
               std::vector<arma::uword>& child_left,
               std::vector<arma::vec>& coef_values,
               std::vector<arma::uvec>& coef_indices,
               std::vector<arma::vec>& leaf_pred_indx,
               std::vector<arma::vec>& leaf_pred_prob,
               std::vector<arma::vec>& leaf_pred_chaz,
               std::vector<double>& leaf_summary,
               arma::vec* pred_horizon);

  void resize_leaves(arma::uword new_size) override;

  double compute_max_leaves() override;

  bool is_col_splittable(arma::uword j) override;

  bool is_node_splittable_internal() override;

  void find_all_cuts() override;

  double compute_split_score() override;

  double compute_mortality(arma::mat& leaf_data);

  void sprout_leaf_internal(uword node_id) override;

  arma::uword predict_value_internal(arma::uvec& pred_leaf_sort,
                                     arma::mat& pred_output,
                                     arma::vec& pred_denom,
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

  void set_unique_event_times(arma::vec event_times){
   this->unique_event_times = &event_times;
  }

  void set_leaf_min_events(double value){
   this->leaf_min_events = value;
  }

  arma::uword find_safe_mtry() override;

  double compute_prediction_accuracy_internal(arma::vec& preds) override;

  arma::mat glm_fit() override;

  // indx holds the times
  std::vector<arma::vec> leaf_pred_indx;
  // prob holds the predicted survival
  std::vector<arma::vec> leaf_pred_prob;
  // chaz holds the cumulative hazard
  std::vector<arma::vec> leaf_pred_chaz;
  // summary (see Tree.h) holds total mortality

  // pointer to event times in forest
  arma::vec* unique_event_times;

  // prediction times
  arma::vec* pred_horizon;

  double leaf_min_events;
  double split_min_events;

 };

 } // namespace aorsf

#endif /* TREESURVIVAL_H_ */
