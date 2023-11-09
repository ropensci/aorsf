
//  ForestClassification.h

#ifndef FORESTCLASSIFICATION_H
#define FORESTCLASSIFICATION_H

#include "Data.h"
#include "globals.h"
#include "Forest.h"

namespace aorsf {

class ForestClassification: public Forest {

public:

 ForestClassification();

 ForestClassification(arma::uword n_class);

 virtual ~ForestClassification() override = default;

 ForestClassification(const ForestClassification&) = delete;
 ForestClassification& operator=(const ForestClassification&) = delete;

 void load(
   arma::uword n_tree,
   arma::uword n_obs,
   std::vector<arma::uvec>& forest_rows_oobag,
   std::vector<std::vector<double>>& forest_cutpoint,
   std::vector<std::vector<arma::uword>>& forest_child_left,
   std::vector<std::vector<arma::vec>>& forest_coef_values,
   std::vector<std::vector<arma::uvec>>& forest_coef_indices,
   std::vector<std::vector<arma::vec>>& forest_leaf_pred_prob,
   std::vector<std::vector<double>>& forest_leaf_summary,
   PartialDepType pd_type,
   std::vector<arma::mat>& pd_x_vals,
   std::vector<arma::uvec>& pd_x_cols,
   arma::vec& pd_probs
 );

 void resize_pred_mat_internal(arma::mat& p) override;

 void compute_prediction_accuracy_internal(
   arma::mat& y,
   arma::vec& w,
   arma::mat& predictions,
   arma::uword row_fill
 ) override;

 // growInternal() in ranger
 void plant() override;

 std::vector<std::vector<arma::vec>> get_leaf_pred_prob();

 uword n_class;

};

}



#endif /* ForestClassification_H */
