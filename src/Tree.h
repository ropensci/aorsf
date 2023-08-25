/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef TREE_H_
#define TREE_H_


#include "Data.h"
#include "globals.h"
#include "utility.h"

 namespace aorsf {

 class Tree {

 public:

  Tree();

  // deleting the copy constructor
  Tree(const Tree&) = delete;
  // deleting the copy assignment operator
  Tree& operator=(const Tree&) = delete;

  void init(Data* data,
            int seed,
            arma::uword mtry,
            double leaf_min_events,
            double leaf_min_obs,
            VariableImportance vi_type,
            double vi_max_pvalue,
            SplitRule split_rule,
            double split_min_events,
            double split_min_obs,
            double split_min_stat,
            arma::uword split_max_cuts,
            arma::uword split_max_retry,
            LinearCombo lincomb_type,
            double lincomb_eps,
            arma::uword lincomb_iter_max,
            bool   lincomb_scale,
            double lincomb_alpha,
            arma::uword lincomb_df_target,
            arma::uword lincomb_ties_method);



  void sample_rows();

  void sample_cols();

  bool is_col_splittable(arma::uword j);

  bool is_node_splittable(arma::uword node_id);

  arma::uvec find_cutpoints();

  double score_logrank();

  double node_split(arma::uvec& cuts_all);

  void node_sprout(arma::uword node_id);

  void grow(arma::vec* vi_numer,
            arma::uvec* vi_denom);

  // void grow(arma::vec& vi_numer, arma::uvec& vi_denom);

  std::vector<arma::uvec>& get_coef_indices() {
   return(coef_indices);
  }

  std::vector<arma::vec>& get_coef_values() {
   return(coef_values);
  }

  std::vector<arma::vec>& get_leaf_pred_horizon(){
   return(leaf_pred_horizon);
  }

  std::vector<arma::vec>& get_leaf_pred_surv(){
   return(leaf_pred_surv);
  }

  std::vector<arma::vec>& get_leaf_pred_chf(){
   return(leaf_pred_chf);
  }

  std::vector<double>& get_cutpoint(){
   return(cutpoint);
  }

  std::vector<arma::uword>& get_child_left(){
   return(child_left);
  }

  // pointers to variable importance in forest
  arma::vec* vi_numer;
  arma::uvec* vi_denom;


  // Pointer to original data
  Data* data;

  arma::uword n_cols_total;
  arma::uword n_rows_total;

  arma::uword n_rows_inbag;

  double n_obs_inbag;
  double n_events_inbag;

  int seed;

  // views of data
  arma::mat x_inbag;
  arma::mat x_oobag;
  arma::mat x_node;

  arma::mat y_inbag;
  arma::mat y_oobag;
  arma::mat y_node;

  // the 'w' is short for 'weights'
  arma::vec w_inbag;
  arma::vec w_oobag;
  arma::vec w_node;

  // g_node indicates where observations will go when this node splits
  // 0 means go down to left node, 1 means go down to right node
  // the 'g' is short for 'groups'
  arma::uvec g_node;

  // which rows of data are held out while growing the tree
  arma::uvec rows_inbag;
  arma::uvec rows_oobag;
  arma::uvec rows_node;
  arma::uvec cols_node;

  // variable importance
  VariableImportance vi_type;
  double vi_max_pvalue;

  // Random number generator
  std::mt19937_64 random_number_generator;

  // tree growing members
  arma::uword mtry;
  double leaf_min_events;
  double leaf_min_obs;

  // node split members
  SplitRule split_rule;
  double split_min_events;
  double split_min_obs;
  double split_min_stat;
  arma::uword split_max_cuts;
  arma::uword split_max_retry;

  // linear combination members
  LinearCombo lincomb_type;
  arma::vec   lincomb;
  arma::uvec  lincomb_sort;
  double      lincomb_eps;
  arma::uword lincomb_iter_max;
  bool        lincomb_scale;
  double      lincomb_alpha;
  arma::uword lincomb_df_target;
  arma::uword lincomb_ties_method;

  // prediction members
  double pred_horizon;

  // predicted values for out-of-bag rows
  arma::vec pred_oobag;

  // which node each inbag observation is currently in.
  arma::uvec node_assignments;

  // cutpoints used to split the nodes
  std::vector<double> cutpoint;

  // left child nodes (right child is left + 1)
  std::vector<arma::uword> child_left;

  // coefficients for linear combinations;
  // one row per variable (mtry rows), one column per node
  // leaf nodes have all coefficients=0
  std::vector<arma::vec> coef_values;
  // std::vector<arma::vec> coef_values;

  // indices of the predictors used by a node
  std::vector<arma::uvec> coef_indices;
  // std::vector<arma::uvec> coef_indices;

  // leaf values (only in leaf nodes)
  std::vector<arma::vec> leaf_pred_horizon;
  std::vector<arma::vec> leaf_pred_surv;
  std::vector<arma::vec> leaf_pred_chf;



 protected:



 };

 } // namespace aorsf

#endif /* TREE_H_ */
