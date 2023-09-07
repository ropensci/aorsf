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

  TreeSurvival(std::vector<double>& cutpoint,
               std::vector<arma::uword>& child_left,
               std::vector<arma::vec>& coef_values,
               std::vector<arma::uvec>& coef_indices,
               std::vector<arma::vec>& leaf_pred_indx,
               std::vector<arma::vec>& leaf_pred_prob,
               std::vector<arma::vec>& leaf_pred_chaz,
               std::vector<double>& leaf_summary);

  double compute_max_leaves() override;

  void resize_leaves(arma::uword new_size) override;

  bool is_col_splittable(arma::uword j) override;

  bool is_node_splittable_internal() override;

  arma::uvec find_cutpoints() override;

  double compute_split_score() override;

  double score_logrank();

  double compute_mortality(arma::mat& leaf_data);

  void node_sprout(uword node_id) override;

  void predict_value(arma::mat* pred_output, arma::vec* pred_denom,
                     arma::vec& pred_times, char pred_type, bool oobag) override;

  std::vector<arma::vec>& get_leaf_pred_indx(){
   return(leaf_pred_indx);
  }

  std::vector<arma::vec>& get_leaf_pred_prob(){
   return(leaf_pred_prob);
  }

  std::vector<arma::vec>& get_leaf_pred_chaz(){
   return(leaf_pred_chaz);
  }

  std::vector<arma::vec> leaf_pred_indx;
  std::vector<arma::vec> leaf_pred_prob;
  std::vector<arma::vec> leaf_pred_chaz;

 };

 } // namespace aorsf

#endif /* TREESURVIVAL_H_ */
