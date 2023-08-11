
//  Forest.h

#ifndef Forest_H
#define Forest_H

#include "Data.h"
#include "globals.h"
#include "Tree.h"

namespace aorsf {

class Forest {

public:

 // Constructor

 Forest();

 // deleting the copy constructor
 Forest(const Forest&) = delete;
 // deleting the copy assignment operator
 Forest& operator=(const Forest&) = delete;

 // Methods

 void init(std::unique_ptr<Data> input_data,
           Rcpp::IntegerVector& tree_seeds,
           arma::uword n_tree,
           arma::uword mtry,
           VariableImportance vi_type,
           // leaves
           double leaf_min_events,
           double leaf_min_obs,
           // node splitting
           SplitRule split_rule,
           double split_min_events,
           double split_min_obs,
           double split_min_stat,
           arma::uword split_max_retry,
           // linear combinations
           LinearCombo lincomb_type,
           double lincomb_eps,
           arma::uword lincomb_iter_max,
           bool lincomb_scale,
           double lincomb_alpha,
           arma::uword lincomb_df_target,
           // predictions
           PredType pred_type,
           bool pred_mode,
           double pred_horizon,
           bool oobag_pred,
           arma::uword oobag_eval_every);

 // virtual void initarma::uwordernal() = 0;

 // Grow or predict
 void run();

 void grow(Function& lincomb_R_function);

 void plant();


 std::vector<std::vector<arma::uvec>> get_coef_indices() {

  std::vector<std::vector<arma::uvec>> result;

  for (auto& tree : trees) {
   result.push_back(tree->get_coef_indices());
  }

  return result;

 }

 // Member variables

 Rcpp::IntegerVector bootstrap_select_times;
 Rcpp::NumericVector bootstrap_select_probs;

 arma::uword n_tree;
 arma::uword mtry;

 Rcpp::IntegerVector tree_seeds;

 std::vector<std::unique_ptr<Tree>> trees;

 std::unique_ptr<Data> data;

 VariableImportance vi_type;

 // leaves
 double leaf_min_events;
 double leaf_min_obs;

 // node splitting
 SplitRule split_rule;
 double split_min_events;
 double split_min_obs;
 double split_min_stat;
 arma::uword split_max_retry;

 // linear combinations
 LinearCombo lincomb_type;
 double      lincomb_eps;
 bool        lincomb_scale;
 double      lincomb_alpha;
 arma::uword lincomb_iter_max;
 arma::uword lincomb_df_target;

 // predictions
 PredType pred_type;
 double   pred_horizon;
 bool     oobag_pred;
 arma::uword oobag_eval_every;


};

}



#endif /* Forest_H */
