
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
           int n_tree,
           int mtry,
           VariableImportance vi_type,
           // leaves
           double leaf_min_events,
           double leaf_min_obs,
           // node splitting
           SplitRule    split_rule,
           double split_min_events,
           double split_min_obs,
           double split_min_stat,
           int    split_max_retry,
           // linear combinations
           LinearCombo  lincomb_type,
           double lincomb_eps,
           int    lincomb_iter_max,
           bool   lincomb_scale,
           double lincomb_alpha,
           int    lincomb_df_target,
           // predictions
           PredType     pred_type,
           double pred_horizon,
           bool   oobag_pred,
           int    oobag_eval_every);

 // virtual void initInternal() = 0;

 // Grow or predict
 void run();

 void grow();

 Rcpp::IntegerVector get_bootstrap_select_times(){
  return bootstrap_select_times;
 }

 Rcpp::NumericVector get_bootstrap_select_probs(){
  return bootstrap_select_probs;
 }


 // Member variables

 Rcpp::IntegerVector bootstrap_select_times;
 Rcpp::NumericVector bootstrap_select_probs;

 int n_tree;
 int mtry;

 Rcpp::IntegerVector tree_seeds;

 std::vector<std::unique_ptr<Tree>> trees;

 std::unique_ptr<Data> data;

 VariableImportance vi_type;

 // leaves
 double leaf_min_events;
 double leaf_min_obs;

 // node splitting
 SplitRule split_rule;
 double    split_min_events;
 double    split_min_obs;
 double    split_min_stat;
 int       split_max_retry;

 // linear combinations
 LinearCombo lincomb_type;
 double      lincomb_eps;
 int         lincomb_iter_max;
 bool        lincomb_scale;
 double      lincomb_alpha;
 int         lincomb_df_target;

 // predictions
 PredType pred_type;
 double   pred_horizon;
 bool     oobag_pred;
 int      oobag_eval_every;


};

}



#endif /* Forest_H */
