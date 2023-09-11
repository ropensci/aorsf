
//  Forest.h

#ifndef FORESTSURVIVAL_H
#define FORESTSURVIVAL_H

#include "Data.h"
#include "globals.h"
#include "Forest.h"

namespace aorsf {

class ForestSurvival: public Forest {

public:

 ForestSurvival();

 ForestSurvival(double leaf_min_events,
                double split_min_events,
                arma::vec& pred_horizon,
                arma::vec& unique_event_times);

 ForestSurvival(const ForestSurvival&) = delete;
 ForestSurvival& operator=(const ForestSurvival&) = delete;

 void load(arma::uword n_tree,
           std::vector<std::vector<double>>& forest_cutpoint,
           std::vector<std::vector<arma::uword>>& forest_child_left,
           std::vector<std::vector<arma::vec>>& forest_coef_values,
           std::vector<std::vector<arma::uvec>>& forest_coef_indices,
           std::vector<std::vector<arma::vec>>& forest_leaf_pred_indx,
           std::vector<std::vector<arma::vec>>& forest_leaf_pred_prob,
           std::vector<std::vector<arma::vec>>& forest_leaf_pred_chaz,
           std::vector<std::vector<double>>& forest_leaf_summary);

 // growInternal() in ranger
 void plant() override;

 std::vector<std::vector<arma::vec>> get_leaf_pred_indx();
 std::vector<std::vector<arma::vec>> get_leaf_pred_prob();
 std::vector<std::vector<arma::vec>> get_leaf_pred_chaz();


};

}



#endif /* Forest_H */
