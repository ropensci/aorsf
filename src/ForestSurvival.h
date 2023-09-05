
//  Forest.h

#ifndef ForestSurvival_H
#define ForestSurvival_H

#include "Data.h"
#include "globals.h"
#include "Forest.h"

namespace aorsf {

class ForestSurvival: public Forest {

public:

 ForestSurvival();

 ForestSurvival(const ForestSurvival&) = delete;
 ForestSurvival& operator=(const ForestSurvival&) = delete;


 void load(arma::uword n_tree,
           std::vector<std::vector<double>>& forest_cutpoint,
           std::vector<std::vector<arma::uword>>& forest_child_left,
           std::vector<std::vector<arma::vec>>& forest_coef_values,
           std::vector<std::vector<arma::uvec>>& forest_coef_indices,
           std::vector<std::vector<arma::vec>>& forest_leaf_pred_horizon,
           std::vector<std::vector<arma::vec>>& forest_leaf_pred_surv,
           std::vector<std::vector<arma::vec>>& forest_leaf_pred_chf,
           std::vector<std::vector<double>>& forest_leaf_pred_mort);


};

}



#endif /* Forest_H */
